var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "AbstractAlgebra.jl",
    "title": "AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#AbstractAlgebra.jl-1",
    "page": "AbstractAlgebra.jl",
    "title": "AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "AbstractAlgebra.jl",
    "title": "Introduction",
    "category": "section",
    "text": "AbstractAlgebra.jl is a computer algebra package for the Julia programming language,  maintained by William Hart, Tommy Hofmann, Claus Fieker and Fredrik Johansson and other interested contributors.Source code\nOnline documentationAbstractAlgebra.jl grew out of the Nemo project after a number of requests from the community for the pure Julia part of Nemo to be split off into a separate project. See the Nemo website for more details about Nemo.Nemo website"
},

{
    "location": "#Features-1",
    "page": "AbstractAlgebra.jl",
    "title": "Features",
    "category": "section",
    "text": "The features of AbstractAlgebra.jl include:Use of Julia multiprecision integers and rationals\nFinite fields (prime order, naive implementation only)\nNumber fields (naive implementation only)\nUnivariate polynomials\nMultivariate polynomials\nRelative and absolute power series\nLaurent series\nFraction fields\nResidue rings, including mathbbZnmathbbZ\nMatrices and linear algebraAll implementations are fully recursive and generic, so that one can build matrices over polynomial rings, over a finite field, for example.AbstractAlgebra.jl also provides a set of abstract types for Groups, Rings, Fields, Modules and elements thereof, which allow external types to be made part of the AbstractAlgebra.jl type hierarchy."
},

{
    "location": "#Installation-1",
    "page": "AbstractAlgebra.jl",
    "title": "Installation",
    "category": "section",
    "text": "To use AbstractAlgebra we require Julia 0.6 or higher. Please see http://julialang.org/downloads for instructions on  how to obtain Julia for your system.At the Julia prompt simply typejulia> Pkg.add(\"AbstractAlgebra\")"
},

{
    "location": "#Quick-start-1",
    "page": "AbstractAlgebra.jl",
    "title": "Quick start",
    "category": "section",
    "text": "Here are some examples of using AbstractAlgebra.jl.This example makes use of multivariate polynomials.using AbstractAlgebra\n\nR, (x, y, z) = PolynomialRing(ZZ, [\"x\", \"y\", \"z\"])\n\nf = x + y + z + 1\n\np = f^20;\n\n@time q = p*(p+1);Here is an example using generic recursive ring constructions.using AbstractAlgebra\n\nR = GF(7)\n\nS, y = PolynomialRing(R, \"y\")\n\nT = ResidueRing(S, y^3 + 3y + 1)\n\nU, z = PolynomialRing(T, \"z\")\n\nf = (3y^2 + y + 2)*z^2 + (2*y^2 + 1)*z + 4y + 3;\n\ng = (7y^2 - y + 7)*z^2 + (3y^2 + 1)*z + 2y + 1;\n\ns = f^4;\n\nt = (s + g)^4;\n\n@time resultant(s, t)Here is an example using matrices.using AbstractAlgebra\n\nR, x = PolynomialRing(ZZ, \"x\")\n\nS = MatrixSpace(R, 10, 10)\n\nM = rand(S, 0:3, -10:10);\n\n@time det(M)And here is an example with power series.using AbstractAlgebra\n\nR, x = QQ[\"x\"]\n\nS, t = PowerSeriesRing(R, 30, \"t\")\n\nu = t + O(t^100)\n\n@time divexact((u*exp(x*u)), (exp(u)-1));"
},

{
    "location": "constructors/#",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing mathematical objects in AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "constructors/#Constructing-mathematical-objects-in-AbstractAlgebra.jl-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing mathematical objects in AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "constructors/#Constructing-objects-in-Julia-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing objects in Julia",
    "category": "section",
    "text": "In Julia, one constructs objects of a given type by calling a type constructor. This is simply a function with the same name as the type itself. For example, to construct a  BigInt object in Julia, we simply call the BigInt constructor:n = BigInt(\"1234567898765434567898765434567876543456787654567890\")Julia also uses constructors to convert between types. For example, to convert an Int to a BigInt:m = BigInt(123)"
},

{
    "location": "constructors/#How-we-construct-objects-in-AbstractAlgebra.jl-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "How we construct objects in AbstractAlgebra.jl",
    "category": "section",
    "text": "As we explain in Appendix A, Julia types don\'t contain enough information to properly model groups, rings, fields, etc. Instead of using types to construct objects, we use special objects that we refer to as parent objects. They behave a lot like Julia types.Consider the following simple example, to create a multiprecision integer:n = ZZ(\"12345678765456787654567890987654567898765678909876567890\")Here ZZ is not a Julia type, but a callable object. However, for most purposes one can think of such a parent object as though it were a type."
},

{
    "location": "constructors/#Constructing-parent-objects-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "Constructing parent objects",
    "category": "section",
    "text": "For more complicated groups, rings, fields, etc., one first needs to construct the parent object before one can use it to construct element objects.AbstractAlgebra.jl provides a set of functions for constructing such parent objects. For example, to create a parent object for univariate polynomials over the integers, we use the PolynomialRing parent object constructor.R, x = PolynomialRing(ZZ, \"x\")\nf = x^3 + 3x + 1\ng = R(12)In this example, R is the parent object and we use it to convert the Int value 12 to an element of the polynomial ring mathbbZx."
},

{
    "location": "constructors/#List-of-parent-object-constructors-1",
    "page": "Constructing mathematical objects in AbstractAlgebra.jl",
    "title": "List of parent object constructors",
    "category": "section",
    "text": "For convenience, we provide a list of all the parent object constructors in AbstractAlgebra.jl and explain what mathematical domains they represent.Mathematics AbstractAlgebra.jl constructor\nR = mathbbZ R = ZZ\nR = mathbbQ R = QQ\nR = mathbbF_p R = GF(p)\nR = mathbbZnmathbbZ R = ResidueRing(ZZ, n)\nS = Rx S, x = PolynomialRing(R, \"x\")\nS = Rx y S, (x, y) = PolynomialRing(R, [\"x\", \"y\"])\nS = Rx (to precision n) S, x = PowerSeriesRing(R, n, \"x\")\nS = R((x)) (to precision n) S, x = LaurentSeriesRing(R, n, \"x\")\nS = K((x)) (to precision n) S, x = LaurentSeriesField(K, n, \"x\")\nS = mboxFrac_R S = FractionField(R)\nS = R(f) S = ResidueRing(R, f)\nS = R(f) (with (f) maximal) S = ResidueField(R, f)\nS = mboxMat_mtimes n(R) S = MatrixSpace(R, m, n)\nS = mathbbQx(f) S, a = NumberField(f, \"a\")"
},

{
    "location": "rings/#",
    "page": "Ring Interface",
    "title": "Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "rings/#Ring-Interface-1",
    "page": "Ring Interface",
    "title": "Ring Interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all rings. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their rings.In addition to the required functions, there are also optional functions which can be provided for certain types of rings, e.g. GCD domains or fields, etc. If implemented, these allow the generic code to provide additional functionality for those rings, or in some cases, to select more efficient algorithms."
},

{
    "location": "rings/#Types-1",
    "page": "Ring Interface",
    "title": "Types",
    "category": "section",
    "text": "Most rings must supply two types:a type for the parent object (representing the ring itself)\na type for elements of that ringFor example, the generic univariate polynomial type in AbstractAlgebra.jl provides two  types in generic/GenericTypes.jl: Generic.PolyRing{T} for the parent objects\nGeneric.Poly{T} for the actual polynomialsThe parent type must belong to AbstractAlgebra.Ring and the element type must belong to AbstractAlgebra.RingElem. Of course, the types may belong to these abstract types transitively, e.g. Poly{T} actually belongs to AbstractAlgebra.PolyElem{T} which in turn belongs to AbstractAlgebra.RingElem.For parameterised rings, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring (see the function base_ring below for a definition).There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "rings/#Parent-object-caches-1",
    "page": "Ring Interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each ring. This means that if the same ring is constructed twice, elements of the two rings will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of rings are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.PolyRing parent objects are looked up in a dictionary PolyID to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn\'t desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two rings  are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two rings if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "rings/#Required-functions-for-all-rings-1",
    "page": "Ring Interface",
    "title": "Required functions for all rings",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for rings in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the ring parent object R and MyElem for the type of the elements of the ring.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for rings that implement that additional functionality, and should be documented as such."
},

{
    "location": "rings/#Data-type-and-parent-object-methods-1",
    "page": "Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "parent_type(::Type{MyElem})Returns the type of the corresponding parent object for the given element type. For example, parent_type(Generic.Poly{T}) will return Generic.PolyRing{T}.elem_type(::Type{MyParent})Returns the type of the elements of the ring whose parent object has the given type. This is the inverse of the parent_type function, i.e. elem_type(Generic.PolyRing{T}) will return Generic.Poly{T}.base_ring(R::MyParent)Given a parent object R, representing a ring, this function returns the parent object of any base ring that parameterises this ring. For example, the base ring of the ring of polynomials over the integers would be the integer ring.If the ring is not parameterised by another ring, this function must return Union{}.Note that there is a distinction between a base ring and other kinds of parameters. For example, in the ring mathbbZnmathbbZ, the modulus n is a parameter, but the only base ring is mathbbZ. We consider the ring mathbbZnmathbbZ to have been constructed from the base ring mathbbZ by taking its quotient by a (principal) ideal.parent(f::MyElem)Return the parent object of the given element, i.e. return the ring to which the given element belongs.This is usually stored in a field parent in each ring element. (If the parent objects have mutable struct types, the internal overhead here is just an additional machine  pointer stored in each element of the ring.)For some element types it isn\'t necessary to append the parent object as a field of every element. This is the case when the parent object can be reconstructed just given the type of the elements. For example, this is the case for the ring of integers and in fact for any ring element type that isn\'t parameterised or generic in any way.isdomain_type(::Type{MyElem})Returns true if every element of the given element type (which may be parameterised or an abstract type) necessarily has a parent that is an integral domain, otherwise if this cannot be guaranteed, the function returns false. For example, if MyElem was the type of elements of generic residue rings of a polynomial ring, the answer to the question would depend on the modulus of the residue  ring. Therefore isdomain_type would have to return false, since we cannot guarantee that we are dealing with elements of an integral domain in general. But if the given element type was for rational integers, the answer would be true, since every rational integer has as parent the ring of rational integers, which is an integral domain.Note that this function depends only on the type of an element and cannot access information about the object itself, or its parent.isexact_type(::Type{MyElem})Returns true if every element of the given type is represented exactly. For example, p-adic numbers, real and complex floating point numbers and power series are not exact, as we can only represent them in general with finite truncations. Similarly polynomials and matrices over inexact element types are themselves inexact.Integers, rationals, finite fields and polynomials and matrices over them are always exact.Note that MyElem may be parameterised or an abstract type, in which case every element of every type represented by MyElem must be exact, otherwise the function must return false.Base.hash(f::MyElem, h::UInt)Return a hash for the object f of type UInt. This is used as a hopefully cheap way to distinguish objects that differ arithmetically. If the object has components, e.g. the coefficients of a polynomial or elements of a matrix, these should be hashed recursively, passing the same parameter h to all levels. Each component should then be xor\'d with h before combining the individual component hashes to give the final hash.The hash functions in AbstractAlgebra.jl usually start from some fixed 64 bit hexadecimal  value that has been picked at random by the library author for that type. That is then truncated to fit a UInt (in case the latter is not 64 bits). This ensures that objects that are the same arithmetically (or that have the same components), but have different types (or structures), are unlikely to hash to the same value.deepcopy_internal(f::MyElem, dict::ObjectIdDict)Return a copy of the given element, recursively copying all components of the object.Obviously the parent, if it is stored in the element, should not be copied. The new element should have precisely the same parent as the old object.For types that cannot self-reference themselves anywhere internally, the dict argument may be ignored.In the case that internal self-references are possible, please consult the Julia documentation on how to implement deepcopy_internal."
},

{
    "location": "rings/#Constructors-1",
    "page": "Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "Outer constructors for most AbstractAlgebra types are provided by overloading the call syntax for parent objects. If R is a parent object for a given ring we provide the following constructors.(R::MyParent)()Return the zero object of the given ring.(R::MyParent)(a::Integer)Coerce the given integer into the given ring.(R::MyParent)(a::MyElem)If a belongs to the given ring, the function returns it (without making a copy). Otherwise an error is thrown.For parameterised rings we also require a function to coerce from the base ring into the parent ring.(R::MyParent{T})(a::T) where T <: AbstractAlgebra.RingElemCoerce a into the ring R if a belongs to the base ring of R."
},

{
    "location": "rings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "zero(R::MyParent)Return the zero element of the given ring.one(R::MyParent)Return the multiplicative identity of the given ring.iszero(f::MyElem)Return true if the given element is the zero element of the ring it belongs to.isone(f::MyElem)Return true if the given element is the multiplicative identity of the ring it belongs to."
},

{
    "location": "rings/#Canonicalisation-1",
    "page": "Ring Interface",
    "title": "Canonicalisation",
    "category": "section",
    "text": "canonical_unit(f::MyElem)When fractions are created with two elements of the given type, it is nice to be able to represent them in some kind of canonical form. This is of course not always possible. But for example, fractions of integers can be canonicalised by first removing any common factors of the numerator and denominator, then making the denominator positive.In AbstractAlgebra.jl, the denominator would be made positive by dividing both the numerator and denominator by the canonical unit of the denominator. For a negative denominator, this would be -1.For elements of a field, canonical_unit simply returns the element itself. In general, canonical_unit of an invertible element should be that element. Finally, if a = ub we should have the identity canonical_unit(a) = canonical_unit(u)*canonical_unit(b).For some rings, it is completely impractical to implement this function, in which case it may return 1 in the given ring. The function must however always exist, and always return an element of the ring."
},

{
    "location": "rings/#String-I/O-1",
    "page": "Ring Interface",
    "title": "String I/O",
    "category": "section",
    "text": "show(io::IO, R::MyParent)This should print (to the given IO object), an English description of the parent ring. If the ring is parameterised, it can call the corresponding show function for any rings it depends on.show(io::IO, f::MyElem)This should print a human readable, textual representation of the object (to the given IO object). It can recursively call the corresponding show functions for any of its components.It may be necessary in some cases to print parentheses around components of f or to print signs of components. For these, the following functions will exist for each component or component type.needs_parentheses(f::MyElem)Should returns true if parentheses are needed around this object when printed, e.g. as a coefficient of a polynomial. As an example, non-constant polynomials would need such parentheses if used as coefficients of another polynomial.displayed_with_minus_in_front(f::MyElem)When printing polynomials, a + sign is usually inserted automatically between terms of the polynomial. However, this is not desirable if the coefficient is negative and that negative sign is already printed when the coefficient is printed.This function must return true if f is printed starting with a negative sign. This suppresses the automatic printing of a + sign by polynomial printing functions that are printing f as a coefficient of a term.Note that if needs_parentheses returns true for f, then displayed_with_minus_in_front should always return false for that f, since an automatic + will need to be printed in front of a coefficient that is printed with parentheses.show_minus_one(::Type{MyElem})When printing polynomials, we prefer to print x rather than 1*x if the degree 1 term has coefficient 1. This can be taken care of without any special support.However, we also prefer to print -x rather than -1*x. This requires special support, since -1 in some rings is not printed as -1 (e.g. -1 in mathbbZ3mathbbZ might be printed as 2). In such rings, show_minus_one should return true.If show_minus_one returns true, polynomial printing functions will not print -x for terms of degree 1 with coefficient -1, but will use the printing function of the given type to print the coefficient in that case."
},

{
    "location": "rings/#Unary-operations-1",
    "page": "Ring Interface",
    "title": "Unary operations",
    "category": "section",
    "text": "-(f::MyElem)Returns -f."
},

{
    "location": "rings/#Binary-operations-1",
    "page": "Ring Interface",
    "title": "Binary operations",
    "category": "section",
    "text": "+(f::MyElem, g::MyElem)\n-(f::MyElem, g::MyElem)\n*(f::MyElem, g::MyElem)Returns f + g, f - g or fg, respectively."
},

{
    "location": "rings/#Comparison-1",
    "page": "Ring Interface",
    "title": "Comparison",
    "category": "section",
    "text": "==(f::MyElem, g::MyElem)Returns true if f and g are arithmetically equal. In the case where the two elements are inexact, the function returns true if they agree to the minimum precision of the two.isequal(f::MyElem, g::MyElem)For exact rings, this should return the same thing as == above. For inexact rings, this returns true only if the two elements are arithmetically equal and have the same precision."
},

{
    "location": "rings/#Powering-1",
    "page": "Ring Interface",
    "title": "Powering",
    "category": "section",
    "text": "^(f::MyElem, e::Int)Return f^e. The function should throw a DomainError() if negative exponents don\'t make sense but are passed to the function."
},

{
    "location": "rings/#Exact-division-1",
    "page": "Ring Interface",
    "title": "Exact division",
    "category": "section",
    "text": "divexact(f::MyElem, g::MyElem)Returns fg, though note that Julia uses / for floating point division. Here we mean exact division in the ring, i.e. return q such that f = gq. A DivideError() should be thrown if g is zero. If no exact quotient exists or an impossible inverse is unavoidably encountered, an error should be thrown."
},

{
    "location": "rings/#Unsafe-operators-1",
    "page": "Ring Interface",
    "title": "Unsafe operators",
    "category": "section",
    "text": "To speed up polynomial and matrix arithmetic, it sometimes makes sense to mutate values in place rather than replace them with a newly created object every time they are modified.For this purpose, certain mutating operators are required. In order to support immutable types (struct in Julia) and systems that don\'t have in-place operators, all unsafe operators must return the (ostensibly) mutated value. Only the returned value is used in computations, so this lifts the requirement that the unsafe operators actually mutate the value.Note the exclamation point is a convention, which indicates that the object may be mutated in-place.To make use of these functions, one must be certain that no other references are held to the object being mutated, otherwise those values will also be changed!The results of deepcopy and all arithmetic operations, including powering and division can be assumed to be new objects without other references being held, as can objects returned from constructors.Note that R(a) where R is the ring a belongs to, does not create a new value. For this case, use deepcopy(a).zero!(f::MyElem)Set the value f to zero in place. Return the mutated value.mul!(c::MyElem, a::MyElem, b::MyElem)Set c to the value ab in place. Return the mutated value. Aliasing is permitted.add!(c::MyElem, a::MyElem, b::MyElem)Set c to the value a + b in place. Return the mutated value. Aliasing is permitted.addeq!(a::MyElem, b::MyElem)Set a to a + b in place. Return the mutated value. Aliasing is permitted."
},

{
    "location": "rings/#Random-generation-1",
    "page": "Ring Interface",
    "title": "Random generation",
    "category": "section",
    "text": "The random functions are only used for test code to generate test data. They therefore don\'t need to provide any guarantees on uniformity, and in fact, test values that are known to be a good source of corner cases can be supplied.rand(R::MyParent, v...)Returns a random element in the given ring of the specified size.There can be as many arguments as is necessary to specify the size of the test example which is being produced."
},

{
    "location": "rings/#Promotion-rules-1",
    "page": "Ring Interface",
    "title": "Promotion rules",
    "category": "section",
    "text": "In order for AbstractAlgebra to be able to automatically coerce up towers of rings, certain promotion rules must be defined. For every ring, one wants to be able to coerce integers into the ring. And for any ring constructed over a base ring, one would like to be able to coerce from the base ring into the ring.The promotion rules look a bit different depending on whether the element type is parameterised or not and whether it is built on a base ring.For ring element types MyElem that are neither parameterised nor built over a base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElemFor ring element types MyType that aren\'t parameterised, but which have a base ring with concrete element type T the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElempromote_rule(::Type{MyElem}, ::Type{T}) = MyElemFor ring element types MyElem{T} that are parameterised by the type of elements of the base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}function promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingEle\nment}\n   promote_rule(T, U) == T ? MyElem{T} : Union{}\nend"
},

{
    "location": "rings/#Required-functionality-for-inexact-rings-1",
    "page": "Ring Interface",
    "title": "Required functionality for inexact rings",
    "category": "section",
    "text": ""
},

{
    "location": "rings/#Approximation-(floating-point-and-ball-arithmetic-only)-1",
    "page": "Ring Interface",
    "title": "Approximation (floating point and ball arithmetic only)",
    "category": "section",
    "text": "isapprox(f::MyElem, g::MyElem; atol::Real=sqrt(eps()))This is used by test code that uses rings involving floating point or ball arithmetic. The function should return true if all components of f and g are equal to within the square root of the Julia epsilon, since numerical noise may make an exact comparison impossible.For parameterised rings over an inexact ring, we also require the following ad hoc approximation functionality.isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemisapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemThese notionally coerce the element of the base ring into the parameterised ring and do a full comparison."
},

{
    "location": "rings/#Optional-functionality-1",
    "page": "Ring Interface",
    "title": "Optional functionality",
    "category": "section",
    "text": "Some functionality is difficult or impossible to implement for all rings in the system. If it is provided, additional functionality or performance may become available. Here is a list of all functions that are considered optional and can\'t be relied on by generic functions in the AbstractAlgebra Ring interface.It may be that no algorithm, or no efficient algorithm is known to implement these functions. As these functions are optional, they do not need to exist. Julia will already inform the user that the function has not been implemented if it is called but doesn\'t exist."
},

{
    "location": "rings/#Optional-basic-manipulation-functionality-1",
    "page": "Ring Interface",
    "title": "Optional basic manipulation functionality",
    "category": "section",
    "text": "isunit(f::MyElem)Return true if the given element is a unit in the ring it belongs to. "
},

{
    "location": "rings/#Optional-binary-ad-hoc-operators-1",
    "page": "Ring Interface",
    "title": "Optional binary ad hoc operators",
    "category": "section",
    "text": "By default, ad hoc operations are handled by AbstractAlgebra.jl if they are not defined explicitly, by coercing both operands into the same ring and then performing the required operation.In some cases, e.g. for matrices, this leads to very inefficient behaviour. In such cases, it is advised to implement some of these operators explicitly.It can occasionally be worth adding a separate set of ad hoc binary operators for the type Int, if this can be done more efficiently than for arbitrary Julia Integer types.+(f::MyElem, c::Integer)\n-(f::MyElem, c::Integer)\n*(f::MyElem, c::Integer)+(c::Integer, f::MyElem)\n-(c::Integer, f::MyElem)\n*(c::Integer, f::MyElem)For parameterised types, it is also sometimes more performant to provide explicit ad hoc operators with elements of the base ring.+(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n-(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n*(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem+(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n-(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n*(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "rings/#Optional-ad-hoc-comparisons-1",
    "page": "Ring Interface",
    "title": "Optional ad hoc comparisons",
    "category": "section",
    "text": "==(f::MyElem, c::Integer)==(c::Integer, f::MyElem)==(f::MyElem{T}, c:T) where T <: AbstractAlgebra.RingElem==(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "rings/#Optional-ad-hoc-exact-division-functions-1",
    "page": "Ring Interface",
    "title": "Optional ad hoc exact division functions",
    "category": "section",
    "text": "divexact(a::MyType{T}, b::T) where T <: AbstractAlgebra.RingElemdivexact(a::MyType, b::Integer)"
},

{
    "location": "rings/#Optional-powering-functions-1",
    "page": "Ring Interface",
    "title": "Optional powering functions",
    "category": "section",
    "text": "^(f::MyElem, e::BigInt)In case f cannot explode in size when powered by a very large integer, and it is practical to do so, one may provide this function to support powering with BigInt exponents (or for external modules, any other big integer type)."
},

{
    "location": "rings/#Optional-unsafe-operators-1",
    "page": "Ring Interface",
    "title": "Optional unsafe operators",
    "category": "section",
    "text": "addmul!(c::MyElem, a::MyElem, b::MyElem, t::MyElem)Set c = c + ab in-place. Return the mutated value. The value t should be a temporary of the same type as a, b and c, which can be used arbitrarily by the implementation to speed up the computation. Aliasing between a, b and c is  permitted."
},

{
    "location": "ncrings/#",
    "page": "Noncommutative ring Interface",
    "title": "Noncommutative ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "ncrings/#Noncommutative-ring-Interface-1",
    "page": "Noncommutative ring Interface",
    "title": "Noncommutative ring Interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all noncommutative rings. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their noncommutative rings."
},

{
    "location": "ncrings/#Types-1",
    "page": "Noncommutative ring Interface",
    "title": "Types",
    "category": "section",
    "text": "Most noncommutative rings must supply two types:a type for the parent object (representing the ring itself)\na type for elements of that ringThe parent type must belong to AbstractAlgebra.NCRing and the element type must belong to AbstractAlgebra.NCRingElem. Of course, the types may belong to these abstract types transitively via an intermediate abstract type.For parameterised rings, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring (see the function base_ring below for a definition).There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "ncrings/#Parent-object-caches-1",
    "page": "Noncommutative ring Interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each ring. This means that if the same ring is constructed twice, elements of the two rings will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of rings are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.PolyRing parent objects are looked up in a dictionary PolyID to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn\'t desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two rings  are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two rings if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "ncrings/#Required-functions-for-all-rings-1",
    "page": "Noncommutative ring Interface",
    "title": "Required functions for all rings",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for noncommutative rings in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the ring parent object R and MyElem for the type of the elements of the ring.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for rings that implement that additional functionality, and should be documented as such."
},

{
    "location": "ncrings/#Data-type-and-parent-object-methods-1",
    "page": "Noncommutative ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "parent_type(::Type{MyElem})Returns the type of the corresponding parent object for the given element type. For example, parent_type(Generic.Poly{T}) will return Generic.PolyRing{T}.elem_type(::Type{MyParent})Returns the type of the elements of the ring whose parent object has the given type. This is the inverse of the parent_type function, i.e. elem_type(Generic.PolyRing{T}) will return Generic.Poly{T}.base_ring(R::MyParent)Given a parent object R, representing a ring, this function returns the parent object of any base ring that parameterises this ring.If the ring is not parameterised by another ring, this function must return Union{}.parent(f::MyElem)Return the parent object of the given element, i.e. return the ring to which the given element belongs.This is usually stored in a field parent in each ring element. (If the parent objects have mutable struct types, the internal overhead here is just an additional machine  pointer stored in each element of the ring.)For some element types it isn\'t necessary to append the parent object as a field of every element. This is the case when the parent object can be reconstructed just given the type of the elements. For example, this is the case for the ring of integers and in fact for any ring element type that isn\'t parameterised or generic in any way.isexact_type(::Type{MyElem})Returns true if every element of the given type is represented exactly.Note that MyElem may be parameterised or an abstract type, in which case every element of every type represented by MyElem must be exact, otherwise the function must return false.Base.hash(f::MyElem, h::UInt)Return a hash for the object f of type UInt. This is used as a hopefully cheap way to distinguish objects that differ arithmetically. If the object has components, e.g. the coefficients of a polynomial or elements of a matrix, these should be hashed recursively, passing the same parameter h to all levels. Each component should then be xor\'d with h before combining the individual component hashes to give the final hash.The hash functions in AbstractAlgebra.jl usually start from some fixed 64 bit hexadecimal  value that has been picked at random by the library author for that type. That is then truncated to fit a UInt (in case the latter is not 64 bits). This ensures that objects that are the same arithmetically (or that have the same components), but have different types (or structures), are unlikely to hash to the same value.deepcopy_internal(f::MyElem, dict::ObjectIdDict)Return a copy of the given element, recursively copying all components of the object.Obviously the parent, if it is stored in the element, should not be copied. The new element should have precisely the same parent as the old object.For types that cannot self-reference themselves anywhere internally, the dict argument may be ignored.In the case that internal self-references are possible, please consult the Julia documentation on how to implement deepcopy_internal."
},

{
    "location": "ncrings/#Constructors-1",
    "page": "Noncommutative ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "Outer constructors for most AbstractAlgebra types are provided by overloading the call syntax for parent objects. If R is a parent object for a given ring we provide the following constructors.(R::MyParent)()Return the zero object of the given ring.(R::MyParent)(a::Integer)Coerce the given integer into the given ring.(R::MyParent)(a::MyElem)If a belongs to the given ring, the function returns it (without making a copy). Otherwise an error is thrown.For parameterised rings we also require a function to coerce from the base ring into the parent ring.(R::MyParent{T})(a::T) where T <: AbstractAlgebra.RingElemCoerce a into the ring R if a belongs to the base ring of R."
},

{
    "location": "ncrings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Noncommutative ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "zero(R::MyParent)Return the zero element of the given ring.one(R::MyParent)Return the multiplicative identity of the given ring.iszero(f::MyElem)Return true if the given element is the zero element of the ring it belongs to.isone(f::MyElem)Return true if the given element is the multiplicative identity of the ring it belongs to."
},

{
    "location": "ncrings/#String-I/O-1",
    "page": "Noncommutative ring Interface",
    "title": "String I/O",
    "category": "section",
    "text": "show(io::IO, R::MyParent)This should print (to the given IO object), an English description of the parent ring. If the ring is parameterised, it can call the corresponding show function for any rings it depends on.show(io::IO, f::MyElem)This should print a human readable, textual representation of the object (to the given IO object). It can recursively call the corresponding show functions for any of its components.It may be necessary in some cases to print parentheses around components of f or to print signs of components. For these, the following functions will exist for each component or component type.needs_parentheses(f::MyElem)Should returns true if parentheses are needed around this object when printed, e.g. as a coefficient of a polynomial. As an example, non-constant polynomials would need such parentheses if used as coefficients of another polynomial.displayed_with_minus_in_front(f::MyElem)When printing polynomials, a + sign is usually inserted automatically between terms of the polynomial. However, this is not desirable if the coefficient is negative and that negative sign is already printed when the coefficient is printed.This function must return true if f is printed starting with a negative sign. This suppresses the automatic printing of a + sign by polynomial printing functions that are printing f as a coefficient of a term.Note that if needs_parentheses returns true for f, then displayed_with_minus_in_front should always return false for that f, since an automatic + will need to be printed in front of a coefficient that is printed with parentheses.show_minus_one(::Type{MyElem})When printing polynomials, we prefer to print x rather than 1*x if the degree 1 term has coefficient 1. This can be taken care of without any special support.However, we also prefer to print -x rather than -1*x. This requires special support, since -1 in some rings is not printed as -1 (e.g. -1 in mathbbZ3mathbbZ might be printed as 2). In such rings, show_minus_one should return true.If show_minus_one returns true, polynomial printing functions will not print -x for terms of degree 1 with coefficient -1, but will use the printing function of the given type to print the coefficient in that case."
},

{
    "location": "ncrings/#Unary-operations-1",
    "page": "Noncommutative ring Interface",
    "title": "Unary operations",
    "category": "section",
    "text": "-(f::MyElem)Returns -f."
},

{
    "location": "ncrings/#Binary-operations-1",
    "page": "Noncommutative ring Interface",
    "title": "Binary operations",
    "category": "section",
    "text": "+(f::MyElem, g::MyElem)\n-(f::MyElem, g::MyElem)\n*(f::MyElem, g::MyElem)Returns f + g, f - g or fg, respectively."
},

{
    "location": "ncrings/#Comparison-1",
    "page": "Noncommutative ring Interface",
    "title": "Comparison",
    "category": "section",
    "text": "==(f::MyElem, g::MyElem)Returns true if f and g are arithmetically equal. In the case where the two elements are inexact, the function returns true if they agree to the minimum precision of the two.isequal(f::MyElem, g::MyElem)For exact rings, this should return the same thing as == above. For inexact rings, this returns true only if the two elements are arithmetically equal and have the same precision."
},

{
    "location": "ncrings/#Powering-1",
    "page": "Noncommutative ring Interface",
    "title": "Powering",
    "category": "section",
    "text": "^(f::MyElem, e::Int)Return f^e. The function should throw a DomainError() if negative exponents don\'t make sense but are passed to the function."
},

{
    "location": "ncrings/#Exact-division-1",
    "page": "Noncommutative ring Interface",
    "title": "Exact division",
    "category": "section",
    "text": "divexact_left(f::MyElem, g::MyElem)\ndivexact_right(f::MyElem, g::MyElem)If f = ga for some a in the ring, the function divexact_left(f, g) returns a. If f = ag then divexact_right(f, g) returns a. A DivideError() should be thrown if division is by zero. If no exact quotient exists or an impossible inverse is unavoidably encountered, an error should be thrown."
},

{
    "location": "ncrings/#Unsafe-operators-1",
    "page": "Noncommutative ring Interface",
    "title": "Unsafe operators",
    "category": "section",
    "text": "To speed up arithmetic of objects built over a given ring, it sometimes makes sense to mutate values in place rather than replace them with a newly created object every time they are modified.For this purpose, certain mutating operators are required. In order to support immutable types (struct in Julia) and systems that don\'t have in-place operators, all unsafe operators must return the (ostensibly) mutated value. Only the returned value is used in computations, so this lifts the requirement that the unsafe operators actually mutate the value.Note the exclamation point is a convention, which indicates that the object may be mutated in-place.To make use of these functions, one must be certain that no other references are held to the object being mutated, otherwise those values will also be changed!The results of deepcopy and all arithmetic operations, including powering and division can be assumed to be new objects without other references being held, as can objects returned from constructors.Note that R(a) where R is the ring a belongs to, does not create a new value. For this case, use deepcopy(a).zero!(f::MyElem)Set the value f to zero in place. Return the mutated value.mul!(c::MyElem, a::MyElem, b::MyElem)Set c to the value ab in place. Return the mutated value. Aliasing is permitted.add!(c::MyElem, a::MyElem, b::MyElem)Set c to the value a + b in place. Return the mutated value. Aliasing is permitted.addeq!(a::MyElem, b::MyElem)Set a to a + b in place. Return the mutated value. Aliasing is permitted."
},

{
    "location": "ncrings/#Random-generation-1",
    "page": "Noncommutative ring Interface",
    "title": "Random generation",
    "category": "section",
    "text": "The random functions are only used for test code to generate test data. They therefore don\'t need to provide any guarantees on uniformity, and in fact, test values that are known to be a good source of corner cases can be supplied.rand(R::MyParent, v...)Returns a random element in the given ring of the specified size.There can be as many arguments as is necessary to specify the size of the test example which is being produced."
},

{
    "location": "ncrings/#Promotion-rules-1",
    "page": "Noncommutative ring Interface",
    "title": "Promotion rules",
    "category": "section",
    "text": "In order for AbstractAlgebra to be able to automatically coerce up towers of rings, certain promotion rules must be defined. For every ring, one wants to be able to coerce integers into the ring. And for any ring constructed over a base ring, one would like to be able to coerce from the base ring into the ring.The promotion rules look a bit different depending on whether the element type is parameterised or not and whether it is built on a base ring.For ring element types MyElem that are neither parameterised nor built over a base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{T}) where {T <: Integer} = MyElemFor ring element types MyType that aren\'t parameterised, but which have a base ring with concrete element type T the promotion rules can be defined as follows:promote_rule(::Type{MyElem}, ::Type{U}) where U <: Integer = MyElempromote_rule(::Type{MyElem}, ::Type{T}) = MyElemFor ring element types MyElem{T} that are parameterised by the type of elements of the base ring, the promotion rules can be defined as follows:promote_rule(::Type{MyElem{T}}, ::Type{MyElem{T}}) where T <: RingElement = MyElem{T}function promote_rule(::Type{MyElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingEle\nment}\n   promote_rule(T, U) == T ? MyElem{T} : Union{}\nend"
},

{
    "location": "ncrings/#Required-functionality-for-inexact-rings-1",
    "page": "Noncommutative ring Interface",
    "title": "Required functionality for inexact rings",
    "category": "section",
    "text": ""
},

{
    "location": "ncrings/#Approximation-(floating-point-and-ball-arithmetic-only)-1",
    "page": "Noncommutative ring Interface",
    "title": "Approximation (floating point and ball arithmetic only)",
    "category": "section",
    "text": "isapprox(f::MyElem, g::MyElem; atol::Real=sqrt(eps()))This is used by test code that uses rings involving floating point or ball arithmetic. The function should return true if all components of f and g are equal to within the square root of the Julia epsilon, since numerical noise may make an exact comparison impossible.For parameterised rings over an inexact ring, we also require the following ad hoc approximation functionality.isapprox(f::MyElem{T}, g::T; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemisapprox(f::T, g::MyElem{T}; atol::Real=sqrt(eps())) where T <: AbstractAlgebra.RingElemThese notionally coerce the element of the base ring into the parameterised ring and do a full comparison."
},

{
    "location": "ncrings/#Optional-functionality-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional functionality",
    "category": "section",
    "text": "Some functionality is difficult or impossible to implement for all rings in the system. If it is provided, additional functionality or performance may become available. Here is a list of all functions that are considered optional and can\'t be relied on by generic functions in the AbstractAlgebra Ring interface.It may be that no algorithm, or no efficient algorithm is known to implement these functions. As these functions are optional, they do not need to exist. Julia will already inform the user that the function has not been implemented if it is called but doesn\'t exist."
},

{
    "location": "ncrings/#Optional-basic-manipulation-functionality-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional basic manipulation functionality",
    "category": "section",
    "text": "isunit(f::MyElem)Return true if the given element is a unit in the ring it belongs to. "
},

{
    "location": "ncrings/#Optional-binary-ad-hoc-operators-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional binary ad hoc operators",
    "category": "section",
    "text": "By default, ad hoc operations are handled by AbstractAlgebra.jl if they are not defined explicitly, by coercing both operands into the same ring and then performing the required operation.In some cases, e.g. for matrices, this leads to very inefficient behaviour. In such cases, it is advised to implement some of these operators explicitly.It can occasionally be worth adding a separate set of ad hoc binary operators for the type Int, if this can be done more efficiently than for arbitrary Julia Integer types.+(f::MyElem, c::Integer)\n-(f::MyElem, c::Integer)\n*(f::MyElem, c::Integer)+(c::Integer, f::MyElem)\n-(c::Integer, f::MyElem)\n*(c::Integer, f::MyElem)For parameterised types, it is also sometimes more performant to provide explicit ad hoc operators with elements of the base ring.+(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n-(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem\n*(f::MyElem{T}, c::T) where T <: AbstractAlgebra.RingElem+(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n-(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem\n*(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "ncrings/#Optional-ad-hoc-comparisons-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional ad hoc comparisons",
    "category": "section",
    "text": "==(f::MyElem, c::Integer)==(c::Integer, f::MyElem)==(f::MyElem{T}, c:T) where T <: AbstractAlgebra.RingElem==(c::T, f::MyElem{T}) where T <: AbstractAlgebra.RingElem"
},

{
    "location": "ncrings/#Optional-powering-functions-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional powering functions",
    "category": "section",
    "text": "^(f::MyElem, e::BigInt)In case f cannot explode in size when powered by a very large integer, and it is practical to do so, one may provide this function to support powering with BigInt exponents (or for external modules, any other big integer type)."
},

{
    "location": "ncrings/#Optional-unsafe-operators-1",
    "page": "Noncommutative ring Interface",
    "title": "Optional unsafe operators",
    "category": "section",
    "text": "addmul!(c::MyElem, a::MyElem, b::MyElem, t::MyElem)Set c = c + ab in-place. Return the mutated value. The value t should be a temporary of the same type as a, b and c, which can be used arbitrarily by the implementation to speed up the computation. Aliasing between a, b and c is  permitted."
},

{
    "location": "euclidean/#",
    "page": "Euclidean Ring Interface",
    "title": "Euclidean Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "euclidean/#Euclidean-Ring-Interface-1",
    "page": "Euclidean Ring Interface",
    "title": "Euclidean Ring Interface",
    "category": "section",
    "text": "If a ring provides a meaningful Euclidean structure such that a useful Euclidean remainder can be computed practically, various additional functionality is provided by AbstractAlgebra.jl for those rings. This functionality depends on the following functions existing.mod(f::MyElem, g::MyElem)Returns the Euclidean remainder of f by g. A DivideError() should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.divrem(f::MyElem, g::MyElem)Returns a pair q, r consisting of the Euclidean quotient and remainder of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.div(f::MyElem, g::MyElem)Returns the Euclidean quotient of f by g. A DivideError should be thrown if g is zero. An error should be thrown if an impossible inverse is encountered.mulmod(f::MyElem, g::MyElem, m::MyElem)Returns fg pmodm.powmod(f::MyElem, e::Int, m::MyElem)Returns f^e pmodm.invmod(f::MyElem, m::MyElem)Returns the inverse of f modulo m. If such an inverse doesn\'t exist, an impossible inverse error should be thrown.divides(f::MyElem, g::MyElem)Returns a pair, flag, q, where flag is set to true if g divides f, in which case the quotient is set to the quotient, or flag is set to false and the quotient is set to zero in the same ring as f and g.remove(f::MyElem, p::MyElem)Returns a pair v, q where p^v is the highest power of p dividing f, and q is the cofactor after f is divided by this power.valuation(f::MyElem, p::MyElem)Returns v where p^v is the highest power of p dividing f.gcd(f::MyElem, g::MyElem)Returns a greatest common divisor of f and g.lcm(f::MyElem, g::MyElem)Returns fggcd(f g) if either f or g is not zero, otherwise it throws a DivideError().gcdx(f::MyElem, g::MyElem)Returns a triple d, s, t such that d = gcd(f g) and d = sf + tg, with s reduced modulo g and t reduced modulo f.gcdinv(f::MyElem, g::MyElem)Returns a tuple d, s such that d = gcd(f g) and s = (fd)^-1 pmodgd. Note that d = 1 iff f is invertible modulo g, in which case s = f^-1 pmodg."
},

{
    "location": "integer/#",
    "page": "Integer ring",
    "title": "Integer ring",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "integer/#Integer-ring-1",
    "page": "Integer ring",
    "title": "Integer ring",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/Integer.jl for making Julia BigInts conform to the AbstractAlgebra.jl Ring interface.In addition to providing a parent object ZZ for Julia BigInts, we implement any additional functionality required by AbstractAlgebra.jl.Because BigInt cannot be directly included in the AbstractAlgebra.jl abstract type hierarchy, we achieve integration of Julia BigInts by introducing a type union, called RingElement, which is a union of AbstractAlgebra.RingElem and a number of Julia types, including BigInt. Everywhere that RingElem is notionally used in AbstractAlgebra.jl, we are in fact using RingElement, with additional care being taken to avoid ambiguities.The details of how this is done are technical, and we refer the reader to the implementation for details. For most intents and purposes, one can think of the Julia BigInt type as belonging to AbstractAlgebra.RingElem.One other technicality is that Julia defines certain functions for BigInt, such as sqrt and exp differently to what AbstractAlgebra.jl requires. To get around this, we redefine these functions internally to AbstractAlgebra.jl, without redefining them for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl to function correctly, without broadcasting pirate definitions of already defined Julia functions to the world.To access the internal definitions, one can use AbstractAlgebra.sqrt and AbstractAlgebra.exp, etc."
},

{
    "location": "integer/#Types-and-parent-objects-1",
    "page": "Integer ring",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Integers have type BigInt, as in Julia itself. We simply supplement the functionality for this type as required for computer algebra.The parent objects of such integers has type Integers{BigInt}.For convenience, we also make Int a part of the AbstractAlgebra.jl type hierarchy and its parent object (accessible as zz) has type Integers{Int}. But we caution that this type is not particularly useful as a model of the integers and may not function as expected within AbstractAlgebra.jl."
},

{
    "location": "integer/#Integer-constructors-1",
    "page": "Integer ring",
    "title": "Integer constructors",
    "category": "section",
    "text": "In order to construct integers in AbstractAlgebra.jl, one can first construct the integer ring itself. This is accomplished using the following constructor.Integers{BigInt}()This gives the unique object of type Integers{BigInt} representing the ring of integers in AbstractAlgebra.jl.In practice, one simply uses ZZ which is assigned to be the return value of the above constructor. There is no need to call the constructor in practice.Here are some examples of creating the integer ring and making use of the resulting parent object to coerce various elements into the ring.Examplesf = ZZ()\ng = ZZ(123)\nh = ZZ(BigInt(1234))"
},

{
    "location": "integer/#Basic-ring-functionality-1",
    "page": "Integer ring",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "The integer ring in AbstractAlgebra.jl implements the full Ring interface and the  Euclidean Ring interface.We give some examples of such functionality.Examplesf = ZZ(12)\n\nh = zero(ZZ)\nk = one(ZZ)\nisone(k) == true\niszero(f) == false\nU = base_ring(ZZ)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)\ng = f + 12\nh = powmod(f, 12, ZZ(17))\nflag, q = divides(f, ZZ(3))"
},

{
    "location": "integer/#Integer-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Integer ring",
    "title": "Integer functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below supplements that provided by Julia itself for its BigInt type."
},

{
    "location": "integer/#AbstractAlgebra.Generic.isunit-Tuple{Integer}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::Integer)\n\nReturn true if a is 1 or -1.\n\n\n\n"
},

{
    "location": "integer/#Basic-functionality-1",
    "page": "Integer ring",
    "title": "Basic functionality",
    "category": "section",
    "text": "isunit(::Integer)Examplesr = ZZ(-1)\n\nisunit(r) == true"
},

{
    "location": "integer/#AbstractAlgebra.sqrt-Tuple{BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.sqrt",
    "category": "method",
    "text": "sqrt(a::T) where T <: Integer\n\nReturn the integer square root of a. If a is not a perfect square an error is thrown.\n\n\n\nsqrt(a::Rational{T}) where T <: Integer\n\nReturn the square root of a if it is the square of a rational, otherwise throw an error.\n\n\n\n"
},

{
    "location": "integer/#AbstractAlgebra.exp-Tuple{BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.exp",
    "category": "method",
    "text": "exp(a::T) where T <: Integer\n\nReturn 1 if a = 0, otherwise throw an exception. This function is not generally of use to the user, but is used internally in AbstractAlgebra.jl.\n\n\n\nexp(a::Rational{T}) where T <: Integer\n\nReturn 1 if a = 0, otherwise throw an exception.\n\n\n\n"
},

{
    "location": "integer/#Square-root-1",
    "page": "Integer ring",
    "title": "Square root",
    "category": "section",
    "text": "AbstractAlgebra.sqrt(a::BigInt)AbstractAlgebra.exp(a::BigInt)Examplesd = AbstractAlgebra.sqrt(ZZ(36))\nm = AbstractAlgebra.exp(ZZ(0))"
},

{
    "location": "integer/#AbstractAlgebra.ppio-Tuple{BigInt,BigInt}",
    "page": "Integer ring",
    "title": "AbstractAlgebra.ppio",
    "category": "method",
    "text": "ppio(a::T, b::T)\n\nSplit a into c*d where c = gcd(a b^infty).\n\n\n\n"
},

{
    "location": "integer/#Coprime-bases-1",
    "page": "Integer ring",
    "title": "Coprime bases",
    "category": "section",
    "text": "ppio(a::BigInt, b::BigInt)Examplesc, n = ppio(ZZ(12), ZZ(26))"
},

{
    "location": "polynomial_rings/#",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Univariate Polynomial Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "polynomial_rings/#Univariate-Polynomial-Ring-Interface-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Univariate Polynomial Ring Interface",
    "category": "section",
    "text": "Univariate polynomial rings are supported in AbstractAlgebra, and in addition to the standard Ring interface, numerous additional functions are required to be present for univariate polynomial rings.Univariate polynomial rings can be built over both commutative and noncommutative rings.Univariate polynomial rings over a field are also Euclidean and therefore such rings must implement the Euclidean interface.Since a sparse distributed multivariate format can generally also handle sparse univariate polynomials, the univariate polynomial interface is designed around the assumption that they are dense. This is not a requirement, but it may be easier to use the multivariate interface for sparse univariate types."
},

{
    "location": "polynomial_rings/#Types-and-parents-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for polynomial rings and their elements over a commutative ring:PolyRing{T} is the abstract type for univariate polynomial ring parent types\nPolyElem{T} is the abstract type for univariate polynomial typesSimilarly there are two abstract types for polynomial rings and their elements over a noncommutative ring:NCPolyRing{T} is the abstract type for univariate polynomial ring parent types\nNCPolyElem{T} is the abstract type for univariate polynomial typesWe have that PolyRing{T} <: AbstractAlgebra.Ring and PolyElem{T} <: AbstractAlgebra.RingElem. Similarly we have that NCPolyRing{T} <: AbstractAlgebra.NCRing and  NCPolyElem{T} <: AbstractAlgebra.NCRingElem.Note that the abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the polynomial ring. For example, in the case of mathbbZx the type T would be the type of an integer, e.g. BigInt.If the parent object for such a ring has type MyZX and polynomials in that ring have type MyZXPoly then one would have:MyZX <: PolyRing{BigInt}\nMyZXPoly <: PolyElem{BigInt}Polynomial rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Polynomial rings should at least be distinguished based on their base (coefficient) ring. But if they have the same base ring and symbol (for their variable/generator), they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "polynomial_rings/#Required-functionality-for-univariate-polynomials-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Required functionality for univariate polynomials",
    "category": "section",
    "text": "In addition to the required functionality for the Ring/NCRing interface (and in the case of polynomials over a field, the Euclidean Ring interface), the Polynomial Ring interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a univariate polynomial ring over R (i.e. S = Rx) with parent object S of type MyPolyRing{T}. We also assume the polynomials in the ring have type MyPoly{T}, where T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem or NCRingElem.We describe the functionality below for polynomials over commutative rings, i.e. with element type belonging to RingElem, however similar constructors should be available for element types belonging to NCRingElem instead, if the coefficient ring is noncommutative."
},

{
    "location": "polynomial_rings/#Constructors-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of coefficients, must be available.(S::MyPolyRing{T})(A::Array{T, 1}) where T <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i].(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: AbstractAlgebra.RingElemCreate the polynomial in the given ring whose degree i coefficient is given by A[i]. The elements of the array are assumed to be able to be coerced into the base ring R.(S::MyPolyRing{T})(A::Array{U, 1}) where T <: AbstractAlgebra.RingElem, U <: IntegerCreate the polynomial in the given ring whose degree i coefficient is given by A[i].It may be desirable to have a additional version of the function that accepts an array of Julia Int values  if this can be done more efficiently.ExamplesS, x = PolynomialRing(QQ, \"x\")\n\nf = S(Rational{BigInt}[2, 3, 1])\ng = S(BigInt[1, 0, 4])\nh = S([4, 7, 2, 9])"
},

{
    "location": "polynomial_rings/#Data-type-and-parent-object-methods-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "var(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn a Symbol representing the variable (generator) of the polynomial ring. Note that this is a Symbol not a String, though its string value will usually be used when printing polynomials.symbols(S::MyPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn the array [s] where s	 is aSymbol` representing the variable of the given polynomial ring. This is provided for uniformity with the multivariate interface, where there is more than one variable, and hence an array of symbols.ExamplesS, x = PolynomialRing(QQ, \"x\")\n\nvsym = var(S)\nV = symbols(S)"
},

{
    "location": "polynomial_rings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "length(f::MyPoly{T}) where T <: AbstractAlgebra.RingElemReturn the length of the given polynomial. The length of the zero polynomial is defined to be 0, otherwise the length is the degree plus 1. The return value should be of type Int.set_length!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemThis function must zero any coefficients beyond the requested length n and then set the length of the polynomial to n. This function does not need to normalise the polynomial and is not useful to the user, but is used extensively by the AbstractAlgebra generic functionality.This function mutates the existing polynomial in-place, but does not return the polynomial.coeff(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the coefficient of the polynomial f of degree n. If n is larger than the degree of the polynomial, it should return zero in the coefficient ring. setcoeff!(f::MyPoly{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElemSet the degree n coefficient of f to a. This mutates the polynomial in-place if possible and returns the mutated polynomial (so that immutable types can also be supported). The function must not assume that the polynomial already has space for n + 1 coefficients. The polynomial must be resized if this is not the case.Note that this function is not required to normalise the polynomial and is not necessarily useful to the user, but is used extensively by the generic functionality in AbstractAlgebra.jl. It is for setting raw coefficients in the representation.normalise(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemGiven a polynomial whose length is currently n, including any leading zero coefficients, return the length of the normalised polynomial (either zero or the length of the polynomial with nonzero leading coefficient). Note that the function does not actually perform the normalisation.fit!(f::MyPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial f internally has space for n coefficients. This function must mutate the function in-place if it is mutable. It does not return the mutated polynomial. Immutable types can still be supported by defining this function to do nothing.Some interfaces for C polynomial types automatically manage the internal allocation of polynomials in every function that can be called on them. Explicit adjustment by the generic code in AbstractAlgebra.jl is not required. In such cases, this function can also be defined to do nothing.ExamplesS, x = PolynomialRing(ZZ, \"x\")\n\nf = x^3 + 3x + 1\ng = S(BigInt[1, 2, 0, 1, 0, 0, 0]);\n\nn = length(f)\nc = coeff(f, 1)\nset_length!(g, normalise(g, 7))\ng = setcoeff!(g, 2, BigInt(11))\nfit!(g, 8)\ng = setcoeff!(g, 7, BigInt(4))\n"
},

{
    "location": "polynomial_rings/#Optional-functionality-for-polynomial-rings-1",
    "page": "Univariate Polynomial Ring Interface",
    "title": "Optional functionality for polynomial rings",
    "category": "section",
    "text": "Sometimes parts of the Euclidean Ring interface can and should be implemented for polynomials over a ring that is not necessarily a field.When divisibility testing can be implemented for a polynomial ring over a field, it  should be possible to implement the following functions from the Euclidean Ring interface:divides\nremove\nvaluationWhen the given polynomial ring is a GCD domain, with an effective GCD algorithm, it may be possible to implement the following functions:gcd\nlcmPolynomial rings can optionally implement any part of the generic univariate polynomial functionality provided by AbstractAlgebra.jl, using the same interface. Obviously additional functionality can also be added to that provided by AbstractAlgebra.jl on an ad hoc basis."
},

{
    "location": "polynomial/#",
    "page": "Generic univariate polynomials",
    "title": "Generic univariate polynomials",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "polynomial/#Generic-univariate-polynomials-1",
    "page": "Generic univariate polynomials",
    "title": "Generic univariate polynomials",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Poly.jl for generic polynomials over any commutative ring belonging to the AbstractAlgebra abstract type hierarchy.As well as implementing the Univariate Polynomial interface, and relevant parts of the Euclidean Ring interface for polynomials over a field, there are many additional generic algorithms implemented for such polynomial rings. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "polynomial/#Types-and-parent-objects-1",
    "page": "Generic univariate polynomials",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Polynomials implemented using the AbstractAlgebra generics have type Generic.Poly{T} where T is the type of elements of the coefficient ring. Internally they consist of a Julia array of coefficients and some additional fields for length and a parent object, etc. See the file src/generic/GenericTypes.jl for details.Parent objects of such polynomials have type Generic.PolyRing{T}.The string representation of the variable of the polynomial ring, and the base/coefficient ring R is stored in the parent object. The polynomial element types belong to the abstract type AbstractAlgebra.PolyElem{T} and the polynomial ring types belong to the abstract type AbstractAlgebra.PolyRing{T}. This enables one to write generic functions that can accept any AbstractAlgebra polynomial type.Note that both the generic polynomial ring type Generic.PolyRing{T} and the abstract type it belongs to, AbstractAlgebra.PolyRing{T} are both called PolyRing. The  former is a (parameterised) concrete type for a polynomial ring over a given base ring whose elements have type T. The latter is an abstract type representing all polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "polynomial/#Polynomial-ring-constructors-1",
    "page": "Generic univariate polynomials",
    "title": "Polynomial ring constructors",
    "category": "section",
    "text": "In order to construct polynomials in AbstractAlgebra.jl, one must first construct the polynomial ring itself. This is accomplished with the following constructor.PolynomialRing(R::AbstractAlgebra.Ring, s::AbstractString; cached::Bool = true)Given a base ring R and string s specifying how the generator (variable) should be printed, return a tuple S, x representing the new polynomial ring S = Rx and the generator x of the ring. By default the parent object S will depend only on R and  x and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.A shorthand version of this function is provided: given a base ring R, we abbreviate the constructor as follows.R[\"x\"]Here are some examples of creating polynomial rings and making use of the resulting parent objects to coerce various elements into the polynomial ring.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = QQ[\"z\"]\n\nf = R()\ng = R(123)\nh = S(BigInt(1234))\nk = S(x + 1)\nm = T(z + 1)All of the examples here are generic polynomial rings, but specialised implementations of polynomial rings provided by external modules will also usually provide a PolynomialRing constructor to allow creation of their polynomial rings."
},

{
    "location": "polynomial/#Basic-ring-functionality-1",
    "page": "Generic univariate polynomials",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "Once a polynomial ring is constructed, there are various ways to construct polynomials in that ring.The easiest way is simply using the generator returned by the PolynomialRing constructor and build up the polynomial using basic arithmetic, as described in the Ring interface. The Julia language also has special syntax for the construction of polynomials in terms of a generator, e.g. we can write 2x instead of 2*x.The polynomial rings in AbstractAlgebra.jl implement the full Ring interface. Of course the entire Univariate Polynomial Ring interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x^3 + 3x + 21\ng = (x + 1)*y^2 + 2x + 1\n\nh = zero(S)\nk = one(R)\nisone(k) == true\niszero(f) == false\nn = length(g)\nU = base_ring(S)\nV = base_ring(y + 1)\nv = var(S)\nT = parent(y + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)For polynomials over a field, the Euclidean Ring interface is implemented.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\nT, y = PolynomialRing(S, \"y\")\n\nf = (3*x^2 + x + 2)*y + x^2 + 1\ng = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1\nh = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1\n\ninvmod(f, g)\nmulmod(f, g, h)\npowmod(f, 3, h)\nh = mod(f, g)\nq, r = divrem(f, g)\nd = gcd(f*h, g*h)\nk = gcdinv(f, h)\nm = lcm(f, h)\nflag, q = divides(g^2, g)\nvaluation(3g^3, g) == 3\nval, q = remove(5g^3, g)\nr, s, t = gcdx(g, h)Functions in the Euclidean Ring interface are supported over residue rings that are not fields, except that if an impossible inverse is encountered during the computation an error is thrown."
},

{
    "location": "polynomial/#Polynomial-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic univariate polynomials",
    "title": "Polynomial functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any polynomial module that implements the full Univariate Polynomial Ring interface. This includes AbstractAlgebra.jl\'s own generic polynomial rings.But if a C library provides all the functionality documented in the Univariate Polynomial Ring interface, then all the functions described here will also be  automatically supplied by AbstractAlgebra.jl for that polynomial type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.modulus-Union{Tuple{PolyElem{T}}, Tuple{T}} where T<:ResElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "method",
    "text": "modulus(a::AbstractAlgebra.PolyElem{T}) where {T <: ResElem}\n\nReturn the modulus of the coefficients of the given polynomial.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.lead-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.lead",
    "category": "method",
    "text": "lead(a::Generic.PolynomialElem)\n\nReturn the leading coefficient of the given polynomial. This will be the nonzero coefficient of the term with highest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.trail-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.trail",
    "category": "method",
    "text": "trail(a::Generic.PolynomialElem)\n\nReturn the trailing coefficient of the given polynomial. This will be the nonzero coefficient of the term with lowest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.gen-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.gen",
    "category": "method",
    "text": "gen(a::MPolyRing{T}, i::Int) where {T <: RingElement}\n\nReturn the i-th generator (variable) of the given polynomial ring.\n\n\n\ngen(R::AbsSeriesRing{T}) where T <: RingElement\n\nReturn the generator of the power series ring, i.e. x + O(x^n) where n is the precision of the power series ring R.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.isgen-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "method",
    "text": "isgen(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is the constant generator of its polynomial ring, otherwise return false.\n\n\n\nisgen(x::MPoly{T}) where {T <: RingElement}\n\nReturn true if the given polynomial is a generator (variable) of the polynomial ring it belongs to.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.isunit-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is a unit in its polynomial ring, otherwise return false.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.ismonomial-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.ismonomial",
    "category": "method",
    "text": "ismonomial(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is a monomial.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.isterm-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.isterm",
    "category": "method",
    "text": "isterm(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is has one term. This function is recursive, with all scalar types returning true.\n\n\n\n"
},

{
    "location": "polynomial/#Basic-functionality-1",
    "page": "Generic univariate polynomials",
    "title": "Basic functionality",
    "category": "section",
    "text": "modulus{T <: ResElem}(::PolyElem{T})lead(::PolyElem)\ntrail(::PolyElem)gen(::PolyElem)isgen(::PolyElem)isunit(::PolyElem)ismonomial(::PolyElem)isterm(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = PolynomialRing(QQ, \"z\")\nU = ResidueRing(ZZ, 17)\nV, w = PolynomialRing(U, \"w\")\n\na = zero(S)\nb = one(S)\n\nc = BigInt(1)//2*z^2 + BigInt(1)//3\nd = x*y^2 + (x + 1)*y + 3\n\nf = lead(d)\ny = gen(S)\ng = isgen(w)\nm = isunit(b)\nn = degree(d)\nr = modulus(w)\nisterm(2y^2) == true\nismonomial(y^2) == true"
},

{
    "location": "polynomial/#Base.truncate-Tuple{PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "Base.truncate",
    "category": "method",
    "text": "truncate(a::Generic.PolynomialElem, n::Int)\n\nReturn a truncated to n terms. \n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.mullow-Union{Tuple{T}, Tuple{PolyElem{T},PolyElem{T},Int64}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.mullow",
    "category": "method",
    "text": "mullow(a::AbstractAlgebra.PolyElem{T}, b::AbstractAlgebra.PolyElem{T}, n::Int) where {T <: RingElement}\n\nReturn atimes b truncated to n terms.\n\n\n\n"
},

{
    "location": "polynomial/#Truncation-1",
    "page": "Generic univariate polynomials",
    "title": "Truncation",
    "category": "section",
    "text": "truncate(::PolyElem, ::Int)mullow{T <: RingElem}(::PolyElem{T}, ::PolyElem{T}, ::Int)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\n\nh = truncate(f, 1)\nk = mullow(f, g, 4)"
},

{
    "location": "polynomial/#Base.reverse-Tuple{PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(x::Generic.PolynomialElem, len::Int)\n\nReturn the reverse of the polynomial x, thought of as a polynomial of the given length (the polynomial will be notionally truncated or padded with zeroes before the leading term if necessary to match the specified length). The resulting polynomial is normalised. If len is negative we throw a DomainError().\n\n\n\n"
},

{
    "location": "polynomial/#Base.reverse-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(x::Generic.PolynomialElem)\n\nReturn the reverse of the polynomial x, i.e. the leading coefficient of x becomes the constant coefficient of the result, etc. The resulting polynomial is normalised.\n\n\n\n"
},

{
    "location": "polynomial/#Reversal-1",
    "page": "Generic univariate polynomials",
    "title": "Reversal",
    "category": "section",
    "text": "reverse(::PolyElem, ::Int)\nreverse(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = reverse(f, 7)\nh = reverse(f)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.shift_left-Tuple{PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.shift_left",
    "category": "method",
    "text": "shift_left(f::Generic.PolynomialElem, n::Int)\n\nReturn the polynomial f shifted left by n terms, i.e. multiplied by x^n.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.shift_right-Tuple{PolyElem,Int64}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.shift_right",
    "category": "method",
    "text": "shift_right(f::Generic.PolynomialElem, n::Int)\n\nReturn the polynomial f shifted right by n terms, i.e. divided by x^n.\n\n\n\n"
},

{
    "location": "polynomial/#Shifting-1",
    "page": "Generic univariate polynomials",
    "title": "Shifting",
    "category": "section",
    "text": "shift_left(::PolyElem, ::Int)shift_right(::PolyElem, ::Int)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = shift_left(f, 7)\nh = shift_right(f, 2)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.pseudorem-Union{Tuple{T}, Tuple{PolyElem{T},PolyElem{T}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.pseudorem",
    "category": "method",
    "text": "pseudorem(f::AbstractAlgebra.PolyElem{T}, g::AbstractAlgebra.PolyElem{T}) where {T <: RingElement}\n\nReturn the pseudoremainder of a divided by b. If b = 0 we throw a DivideError().\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.pseudodivrem-Union{Tuple{T}, Tuple{PolyElem{T},PolyElem{T}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.pseudodivrem",
    "category": "method",
    "text": "pseudodivrem(f::AbstractAlgebra.PolyElem{T}, g::AbstractAlgebra.PolyElem{T}) where {T <: RingElement}\n\nReturn a tuple (q r) consisting of the pseudoquotient and pseudoremainder of a divided by b. If b = 0 we throw a DivideError().\n\n\n\n"
},

{
    "location": "polynomial/#Pseudodivision-1",
    "page": "Generic univariate polynomials",
    "title": "Pseudodivision",
    "category": "section",
    "text": "Given two polynomials a b, pseudodivision computes polynomials q and r with length(r)  length(b) such that L^d a = bq + r where d = length(a) - length(b) + 1 and L is the leading coefficient of b.We call q the pseudoquotient and r the pseudoremainder.pseudorem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})pseudodivrem{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\n\nh = pseudorem(f, g)\nq, r = pseudodivrem(f, g)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.content-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.content",
    "category": "method",
    "text": "content(a::AbstractAlgebra.PolyElem)\n\nReturn the content of a, i.e. the greatest common divisor of its coefficients.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.primpart-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.primpart",
    "category": "method",
    "text": "primpart(a::AbstractAlgebra.PolyElem)\n\nReturn the primitive part of a, i.e. the polynomial divided by its content.\n\n\n\n"
},

{
    "location": "polynomial/#Content-and-primitive-part-1",
    "page": "Generic univariate polynomials",
    "title": "Content and primitive part",
    "category": "section",
    "text": "content(::PolyElem)primpart(::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nk = x*y^2 + (x + 1)*y + 3\n\nn = content(k)\np = primpart(k*(x^2 + 1))"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{T}, Tuple{PolyElem{T},T}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.PolyElem, b::T) where {T <: RingElement}\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.evaluate-Tuple{PolyElem,Integer}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.PolyElem, b::T) where {T <: RingElement}\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.compose-Tuple{PolyElem,PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.compose",
    "category": "method",
    "text": "compose(a::AbstractAlgebra.PolyElem, b::AbstractAlgebra.PolyElem)\n\nCompose the polynomial a with the polynomial b and return the result, i.e. return acirc b.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.subst-Union{Tuple{T}, Tuple{PolyElem{T},Any}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.subst",
    "category": "method",
    "text": "subst(f::AbstractAlgebra.PolyElem{T}, a::Any) where {T <: RingElement}\n\nEvaluate the polynomial f at a. Note that a can be anything, whether a ring element or not.\n\n\n\n"
},

{
    "location": "polynomial/#Evaluation,-composition-and-substitution-1",
    "page": "Generic univariate polynomials",
    "title": "Evaluation, composition and substitution",
    "category": "section",
    "text": "evaluate{T <: RingElem}(::PolyElem{T}, ::T)\nevaluate(::PolyElem, ::Integer)compose(::PolyElem, ::PolyElem)subst{T <: RingElem}(::PolyElem{T}, ::Any)We also overload the functional notation so that the polynomial f can be evaluated at a by writing f(a). ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n   \nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\nM = R[x + 1 2x; x - 3 2x - 1]\n\nk = evaluate(f, 3)\nm = evaluate(f, x^2 + 2x + 1)\nn = compose(f, g)\np = subst(f, M)\nq = f(M)\nr = f(23)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.derivative-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.derivative",
    "category": "method",
    "text": "derivative(a::Generic.PolynomialElem)\n\nReturn the derivative of the polynomial a.\n\n\n\nderivative(f::MPoly{T}, x::MPoly{T}) where {T <: RingElement}\n\nReturn the partial derivative of f with respect to x. The value x must be a generator of the polynomial ring of f.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.integral-Union{Tuple{PolyElem{T}}, Tuple{T}} where T<:Union{FieldElem, ResElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.integral",
    "category": "method",
    "text": "integral(x::AbstractAlgebra.PolyElem{T}) where {T <: Union{AbstractAlgebra.ResElem, FieldElement}}\n\nReturn the integral of the polynomial x.\n\n\n\n"
},

{
    "location": "polynomial/#Derivative-and-integral-1",
    "page": "Generic univariate polynomials",
    "title": "Derivative and integral",
    "category": "section",
    "text": "derivative(::PolyElem)integral{T <: Union{ResElem, FieldElem}}(::PolyElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\nT, z = PolynomialRing(QQ, \"z\")\nU = ResidueRing(T, z^3 + 3z + 1)\nV, w = PolynomialRing(U, \"w\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (z^2 + 2z + 1)*w^2 + (z + 1)*w - 2z + 4\n\nh = derivative(f)\nk = integral(g)   "
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.resultant-Union{Tuple{T}, Tuple{PolyElem{T},PolyElem{T}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.resultant",
    "category": "method",
    "text": "resultant(p::AbstractAlgebra.PolyElem{T}, q::AbstractAlgebra.PolyElem{T}) where {T <: RingElement}\n\nReturn the resultant of the given polynomials.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.resx-Union{Tuple{T}, Tuple{PolyElem{T},PolyElem{T}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.resx",
    "category": "method",
    "text": "resx(a::AbstractAlgebra.PolyElem{T}, b::AbstractAlgebra.PolyElem{T}) where {T <: RingElement}\n\nReturn a tuple (r s t) such that r is the resultant of a and b and such that r = atimes s + btimes t.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.discriminant-Tuple{PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.discriminant",
    "category": "method",
    "text": "discriminant(a::AbstractAlgebra.PolyElem)\n\nReturn the discrimnant of the given polynomial.\n\n\n\n"
},

{
    "location": "polynomial/#Resultant-and-discriminant-1",
    "page": "Generic univariate polynomials",
    "title": "Resultant and discriminant",
    "category": "section",
    "text": "resultant{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})resx{T <: RingElem}(::PolyElem{T}, ::PolyElem{T})discriminant(a::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = 3x*y^2 + (x + 1)*y + 3\ng = 6(x + 1)*y + (x^3 + 2x + 2)\n\nh = resultant(f, g)\nk = discriminant(f)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.monomial_to_newton!-Union{Tuple{T}, Tuple{Array{T,1},Array{T,1}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.monomial_to_newton!",
    "category": "method",
    "text": "monomial_to_newton!(P::Array{T, 1}, roots::Array{T, 1}) where {T <: RingElement}\n\nConverts a polynomial p, given as an array of coefficients, in-place from its coefficients given in the standard monomial basis to the Newton basis for the roots r_0 r_1 ldots r_n-2. In other words, this determines output coefficients c_i such that c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + ldots + c_n-1(x-r_0)(x-r_1)cdots(x-r_n-2) is equal to the input polynomial.\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.newton_to_monomial!-Union{Tuple{T}, Tuple{Array{T,1},Array{T,1}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.newton_to_monomial!",
    "category": "method",
    "text": "newton_to_monomial!(P::Array{T, 1}, roots::Array{T, 1}) where {T <: RingElement}\n\nConverts a polynomial p, given as an array of coefficients, in-place from its coefficients given in the Newton basis for the roots r_0 r_1 ldots r_n-2 to the standard monomial basis. In other words, this evaluates c_0 + c_1(x-r_0) + c_2(x-r_0)(x-r_1) + ldots + c_n-1(x-r_0)(x-r_1)cdots(x-r_n-2) where c_i are the input coefficients given by p.\n\n\n\n"
},

{
    "location": "polynomial/#Newton-representation-1",
    "page": "Generic univariate polynomials",
    "title": "Newton representation",
    "category": "section",
    "text": "monomial_to_newton!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})newton_to_monomial!{T <: RingElem}(::Array{T, 1}, ::Array{T, 1})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = 3x*y^2 + (x + 1)*y + 3\ng = deepcopy(f)\nroots = [R(1), R(2), R(3)]\n\nmonomial_to_newton!(g.coeffs, roots)\nnewton_to_monomial!(g.coeffs, roots)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.interpolate-Union{Tuple{T}, Tuple{PolyRing,Array{T,1},Array{T,1}}} where T<:RingElem",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.interpolate",
    "category": "method",
    "text": "interpolate(S::AbstractAlgebra.PolyRing, x::Array{T, 1}, y::Array{T, 1}) where {T <: RingElement}\n\nGiven two arrays of values xs and ys of the same length n, find the polynomial f in the polynomial ring R of length at most n such that f has the value ys at the points xs. The values in the arrays xs and ys must belong to the base ring of the polynomial ring R. If no such polynomial exists, an exception is raised.\n\n\n\n"
},

{
    "location": "polynomial/#Interpolation-1",
    "page": "Generic univariate polynomials",
    "title": "Interpolation",
    "category": "section",
    "text": "interpolate{T <: RingElem}(::PolyRing, ::Array{T, 1}, ::Array{T, 1})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nxs = [R(1), R(2), R(3), R(4)]\nys = [R(1), R(4), R(9), R(16)]\n\nf = interpolate(S, xs, ys)"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.chebyshev_t-Tuple{Int64,PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.chebyshev_t",
    "category": "method",
    "text": "chebyshev_t(n::Int, x::AbstractAlgebra.PolyElem)\n\nReturn the Chebyshev polynomial of the first kind T_n(x), defined by T_n(x) = cos(n cos^-1(x)).\n\n\n\n"
},

{
    "location": "polynomial/#AbstractAlgebra.Generic.chebyshev_u-Tuple{Int64,PolyElem}",
    "page": "Generic univariate polynomials",
    "title": "AbstractAlgebra.Generic.chebyshev_u",
    "category": "method",
    "text": "chebyshev_u(n::Int, x::AbstractAlgebra.PolyElem)\n\nReturn the Chebyshev polynomial of the first kind U_n(x), defined by (n+1) U_n(x) = T_n+1(x).\n\n\n\n"
},

{
    "location": "polynomial/#Special-functions-1",
    "page": "Generic univariate polynomials",
    "title": "Special functions",
    "category": "section",
    "text": "The following special functions can be computed for any polynomial ring. Typically one uses the generator x of a polynomial ring to get the respective special polynomials expressed in terms of that generator.chebyshev_t(::Int, ::PolyElem)chebyshev_u(::Int, ::PolyElem)ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS, y = PolynomialRing(R, \"y\")\n\nf = chebyshev_t(20, y)\ng = chebyshev_u(15, y)"
},

{
    "location": "ncpolynomial/#",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Generic univariate polynomials over a noncommutative ring",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "ncpolynomial/#Generic-univariate-polynomials-over-a-noncommutative-ring-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Generic univariate polynomials over a noncommutative ring",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/NCPoly.jl for generic polynomials over any noncommutative ring belonging to the AbstractAlgebra abstract type hierarchy.As well as implementing the Univariate Polynomial interface, there are many additional generic algorithms implemented for such polynomial rings. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "ncpolynomial/#Types-and-parent-objects-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Polynomials implemented using the AbstractAlgebra generics have type Generic.NCPoly{T} where T is the type of elements of the coefficient ring. Internally they consist of a Julia array of coefficients and some additional fields for length and a parent object, etc. See the file src/generic/GenericTypes.jl for details.Parent objects of such polynomials have type Generic.NCPolyRing{T}.The string representation of the variable of the polynomial ring, and the base/coefficient ring R is stored in the parent object. The polynomial element types belong to the abstract type AbstractAlgebra.NCPolyElem{T} and the polynomial ring types belong to the abstract type AbstractAlgebra.NCPolyRing{T}. This enables one to write generic functions that can accept any AbstractAlgebra polynomial type.Note that both the generic polynomial ring type Generic.NCPolyRing{T} and the abstract type it belongs to, AbstractAlgebra.NCPolyRing{T} are both called NCPolyRing. The  former is a (parameterised) concrete type for a polynomial ring over a given base ring whose elements have type T. The latter is an abstract type representing all polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "ncpolynomial/#Polynomial-ring-constructors-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Polynomial ring constructors",
    "category": "section",
    "text": "In order to construct polynomials in AbstractAlgebra.jl, one must first construct the polynomial ring itself. This is accomplished with the following constructor.PolynomialRing(R::AbstractAlgebra.NCRing, s::AbstractString; cached::Bool = true)Given a base ring R and string s specifying how the generator (variable) should be printed, return a tuple S, x representing the new polynomial ring S = Rx and the generator x of the ring. By default the parent object S will depend only on R and  x and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.A shorthand version of this function is provided: given a base ring R, we abbreviate the constructor as follows.R[\"x\"]Here are some examples of creating polynomial rings and making use of the resulting parent objects to coerce various elements into the polynomial ring.ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(T, \"y\")\nU, z = R[\"z\"]\n\nf = S()\ng = S(123)\nh = T(BigInt(1234))\nk = T(x + 1)\nm = U(z + 1)All of the examples here are generic polynomial rings, but specialised implementations of polynomial rings provided by external modules will also usually provide a PolynomialRing constructor to allow creation of their polynomial rings."
},

{
    "location": "ncpolynomial/#Basic-ring-functionality-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "Once a polynomial ring is constructed, there are various ways to construct polynomials in that ring.The easiest way is simply using the generator returned by the PolynomialRing constructor and build up the polynomial using basic arithmetic, as described in the Ring interface. The Julia language also has special syntax for the construction of polynomials in terms of a generator, e.g. we can write 2x instead of 2*x.The polynomial rings in AbstractAlgebra.jl implement the full Ring interface. Of course the entire Univariate Polynomial Ring interface is also implemented.We give some examples of such functionality.ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(T, \"y\")\n\nf = x^3 + 3x + 21\ng = (x + 1)*y^2 + 2x + 1\n\nh = zero(T)\nk = one(S)\nisone(k) == true\niszero(f) == false\nn = length(g)\nU = base_ring(T)\nV = base_ring(y + 1)\nv = var(T)\nU = parent(y + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)"
},

{
    "location": "ncpolynomial/#Polynomial-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Polynomial functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any polynomial module that implements the full Univariate Polynomial Ring interface over a noncommutative ring. This includes AbstractAlgebra.jl\'s own generic polynomial rings.But if a C library provides all the functionality documented in the Univariate Polynomial Ring interface over a noncommutative ring, then all the functions described here will also be automatically supplied by AbstractAlgebra.jl for that polynomial type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.lead-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.lead",
    "category": "method",
    "text": "lead(a::Generic.PolynomialElem)\n\nReturn the leading coefficient of the given polynomial. This will be the nonzero coefficient of the term with highest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.trail-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.trail",
    "category": "method",
    "text": "trail(a::Generic.PolynomialElem)\n\nReturn the trailing coefficient of the given polynomial. This will be the nonzero coefficient of the term with lowest degree unless the polynomial in the zero polynomial, in which case a zero coefficient is returned.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.gen-Tuple{AbstractAlgebra.NCPolyRing}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.gen",
    "category": "method",
    "text": "gen(R::AbstractAlgebra.NCPolyRing)\n\nReturn the generator of the given polynomial ring.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.isgen-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "method",
    "text": "isgen(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is the constant generator of its polynomial ring, otherwise return false.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.isunit-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is a unit in its polynomial ring, otherwise return false.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.ismonomial-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.ismonomial",
    "category": "method",
    "text": "ismonomial(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is a monomial.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.isterm-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.isterm",
    "category": "method",
    "text": "isterm(a::Generic.PolynomialElem)\n\nReturn true if the given polynomial is has one term. This function is recursive, with all scalar types returning true.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Basic-functionality-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Basic functionality",
    "category": "section",
    "text": "lead(::NCPolyElem)\ntrail(::NCPolyElem)gen(::NCPolyRing)isgen(::NCPolyElem)isunit(::NCPolyElem)ismonomial(::NCPolyElem)isterm(::NCPolyElem)ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n\na = zero(T)\nb = one(T)\n\nc = BigInt(1)*z^2 + BigInt(1)\nd = x*y^2 + (x + 1)*y + 3\n\nf = lead(d)\ny = gen(T)\ng = isgen(y)\nm = isunit(b)\nn = degree(d)\nisterm(2y^2) == true\nismonomial(y^2) == true"
},

{
    "location": "ncpolynomial/#Base.truncate-Tuple{NCPolyElem,Int64}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Base.truncate",
    "category": "method",
    "text": "truncate(a::Generic.PolynomialElem, n::Int)\n\nReturn a truncated to n terms. \n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.mullow-Union{Tuple{T}, Tuple{NCPolyElem{T},NCPolyElem{T},Int64}} where T<:NCRingElem",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.mullow",
    "category": "method",
    "text": "mullow(a::AbstractAlgebra.NCPolyElem{T}, b::AbstractAlgebra.NCPolyElem{T}, n::Int) where {T <: NCRingElem}\n\nReturn atimes b truncated to n terms.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Truncation-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Truncation",
    "category": "section",
    "text": "truncate(::NCPolyElem, ::Int)mullow(::NCPolyElem{T}, ::NCPolyElem{T}, ::Int) where T <: NCRingElemExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\ng = (x + 1)*y + (x^3 + 2x + 2)\n\nh = truncate(f, 1)\nk = mullow(f, g, 4)"
},

{
    "location": "ncpolynomial/#Base.reverse-Tuple{NCPolyElem,Int64}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(x::Generic.PolynomialElem, len::Int)\n\nReturn the reverse of the polynomial x, thought of as a polynomial of the given length (the polynomial will be notionally truncated or padded with zeroes before the leading term if necessary to match the specified length). The resulting polynomial is normalised. If len is negative we throw a DomainError().\n\n\n\n"
},

{
    "location": "ncpolynomial/#Base.reverse-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Base.reverse",
    "category": "method",
    "text": "reverse(x::Generic.PolynomialElem)\n\nReturn the reverse of the polynomial x, i.e. the leading coefficient of x becomes the constant coefficient of the result, etc. The resulting polynomial is normalised.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Reversal-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Reversal",
    "category": "section",
    "text": "reverse(::NCPolyElem, ::Int)\nreverse(::NCPolyElem)ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = reverse(f, 7)\nh = reverse(f)"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.shift_left-Tuple{NCPolyElem,Int64}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.shift_left",
    "category": "method",
    "text": "shift_left(f::Generic.PolynomialElem, n::Int)\n\nReturn the polynomial f shifted left by n terms, i.e. multiplied by x^n.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.shift_right-Tuple{NCPolyElem,Int64}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.shift_right",
    "category": "method",
    "text": "shift_right(f::Generic.PolynomialElem, n::Int)\n\nReturn the polynomial f shifted right by n terms, i.e. divided by x^n.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Shifting-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Shifting",
    "category": "section",
    "text": "shift_left(::NCPolyElem, ::Int)shift_right(::NCPolyElem, ::Int)ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\ng = shift_left(f, 7)\nh = shift_right(f, 2)"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{T}, Tuple{NCPolyElem{T},T}} where T<:NCRingElem",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.NCPolyElem, b::T) where {T <: NCRingElem}\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.evaluate-Tuple{NCPolyElem,Integer}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.NCPolyElem, b::Union{Integer, Rational, AbstractFloat})\n\nEvaluate the polynomial a at the value b and return the result.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Evaluation-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate{T <: NCRingElem}(::NCPolyElem{T}, ::T)\nevaluate(::NCPolyElem, ::Integer)We also overload the functional notation so that the polynomial f can be evaluated at a by writing f(a). ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n   \nf = x*y^2 + (x + 1)*y + 3\n\nk = evaluate(f, 3)\nm = evaluate(f, x^2 + 2x + 1)\nr = f(23)"
},

{
    "location": "ncpolynomial/#AbstractAlgebra.Generic.derivative-Tuple{NCPolyElem}",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "AbstractAlgebra.Generic.derivative",
    "category": "method",
    "text": "derivative(a::Generic.PolynomialElem)\n\nReturn the derivative of the polynomial a.\n\n\n\n"
},

{
    "location": "ncpolynomial/#Derivative-1",
    "page": "Generic univariate polynomials over a noncommutative ring",
    "title": "Derivative",
    "category": "section",
    "text": "derivative(::NCPolyElem)ExamplesR = MatrixAlgebra(ZZ, 2)\nS, x = PolynomialRing(R, \"x\")\nT, y = PolynomialRing(S, \"y\")\n\nf = x*y^2 + (x + 1)*y + 3\n\nh = derivative(f)"
},

{
    "location": "mpolynomial_rings/#",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Multivariate Polynomial Ring Interface",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "mpolynomial_rings/#Multivariate-Polynomial-Ring-Interface-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Multivariate Polynomial Ring Interface",
    "category": "section",
    "text": "Multivariate polynomial rings are supported in AbstractAlgebra.jl, and in addition to the standard Ring interface, numerous additional functions are provided.Unlike other kinds of rings, even complex operations such as GCD depend heavily on the multivariate representation. Therefore AbstractAlgebra.jl cannot provide much in the way of additional functionality to external multivariate implementations.This means that external libraries must be able to implement their multivariate formats in whatever way they see fit. The required interface here should be implemented, even if it is not optimal. But it can be extended, either by implementing one of the optional interfaces, or by extending the required interface in some other way.Naturally, any multivariate polynomial ring implementation provides the full Ring interface, in order to be treated as a ring for the sake of AbstractAlgebra.jl.Considerations which make it impossible for AbstractAlgebra.jl to provide generic functionality on top of an arbitrary multivariate module include:orderings (lexical, degree, weighted, block, arbitrary)\nsparse or dense representation\ndistributed or recursive representation\npacked or unpacked exponents\nexponent bounds (and whether adaptive or not)\nrandom access or iterators\nwhether monomials and polynomials have the same type\nwhether special cache aware data structures such as Geobuckets are used"
},

{
    "location": "mpolynomial_rings/#Types-and-parents-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra.jl provides two abstract types for multivariate polynomial rings and their elements:MPolyRing{T} is the abstract type for multivariate polynomial ring parent types\nMPolyElem{T} is the abstract type for multivariate polynomial typesWe have that MPolyRing{T} <: AbstractAlgebra.Ring and  MPolyElem{T} <: AbstractAlgebra.RingElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the polynomial ring. For example, in the case of mathbbZx y the type T would be the type of an integer, e.g. BigInt.Multivariate polynomial rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Multivariate polynomial rings should at least be distinguished based on their base (coefficient) ring and number of variables. But if they have the same base ring, symbols (for their variables/generators) and ordering, they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "mpolynomial_rings/#Required-functionality-for-multivariate-polynomials-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Required functionality for multivariate polynomials",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface, the Multivariate Polynomial interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a multivariate polynomial ring over R (i.e. S = Rx y ldots) with parent object S of type MyMPolyRing{T}. We also assume the polynomials in the ring have type MyMPoly{T}, where T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem or more generally the union type RingElement which includes the Julia integer, rational and floating point types."
},

{
    "location": "mpolynomial_rings/#Data-type-and-parent-object-methods-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "symbols(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn an array of Symbols representing the variables (generators) of the polynomial ring. Note that these are Symbols not Strings, though their string values will usually be used when printing polynomials.nvars(f::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn the number of variables of the polynomial ring.gens(S::MyMPolyRing{T}) where T <: AbstractAlgebra.RingElemReturn an array of all the generators (variables) of the given polynomial ring (as polynomials).The first entry in the array will be the variable with most significance with respect to the ordering.gen(S::MyMPolyRing{T}, i::Int) where T <: AbstractAlgebra.RingElemReturn the i-th generator (variable) of the given polynomial ring (as a polynomial).ordering(S::MyMPolyRing{T})Return the ordering of the given polynomial ring as a symbol. Supported values currently include :lex, :deglex and :degrevlex.ExamplesS, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"]; ordering=:deglex)\n\nV = symbols(S)\nX = gens(S)\nord = ordering(S)"
},

{
    "location": "mpolynomial_rings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "length(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the number of nonzero terms of the given polynomial. The length of the zero polynomial is defined to be 0. The return value should be of type Int.degree(f::MyMPoly{T}, i::Int) where T <: AbstractAlgebra.RingElemReturn the degree of the polynomial f in the i-th variable.degrees(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn an array of the degrees of the polynomial f in each of the variables.total_degree(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the total degree of the polynomial f, i.e. the highest sum of exponents occuring in any term of f.isgen(x::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if x is a generator of the polynomial ring.isunit(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f is a unit in its parent polynomial ring.isconstant(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f is a constant polynomial. The zero polynomial is considered constant for the purposes of this function.isterm(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f consists of a single term.ismonomial(f::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn true if f consists of a single term with coefficient 1.coeffs(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn an iterator for the coefficients of the polynomial p, starting with the coefficient of the most significant term with respect to the ordering. Generic code will provide this function automatically for random access polynomials that implement the coeff function.monomials(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn an iterator for the monomials of the polynomial p, starting with the monomial of the most significant term with respect to the ordering. Monomials in AbstractAlgebra are defined to have coefficient 1. See the function terms if you also require the coefficients, however note that only monomials can be compared. Generic code will provide this function automatically for random access polynomials that implement the monomial function.terms(p::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn an iterator for the terms of the polynomial p, starting with the most significant term with respect to the ordering. Terms in AbstractAlgebra include the coefficient. Generic code will provide this function automatically for random access polynomials that implement the term function.exponent_vectors(a::MyMPoly{T}) where T <: RingElementReturn an iterator for the exponent vectors for each of the terms of the polynomial starting with the most significant term with respect to the ordering. Each exponent vector is an array of Ints, one for each variable, in the order given when the polynomial ring was created. Generic code will provide this function automatically for random access polynomials that implement the exponent_vector function.ExamplesS, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^3*y + 3x*y^2 + 1\n\nn = length(f)\nisgen(y) == true\nB, b = max_degrees(f)\nnvars(f) == 2\nisunit(f) == false\nisconstant(f) == false\nisterm(2x*y) == true\nismonomial(x*y) == false\nV = vars(f)\nC = collect(coeffs(f))\nM = collect(monomials(f))\nT = collect(terms(f))\nV = collect(exponent_vectors(f))\ndegree(f, 2) == 2\nd = total_degree(f)"
},

{
    "location": "mpolynomial_rings/#Exact-division-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Exact division",
    "category": "section",
    "text": "For any ring that implements exact division, the following can be implemented.divexact(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn the exact quotient of f by g if it exists, otherwise throw an error.divides(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a tuple (flag, q) where flag is true if g divides f, in which case q will be the exact quotient, or flag is false and q is set to zero.remove(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturns a tuple (v q) such that the highest power of g that divides f is g^v and the cofactor is q.valuation(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturns v such that the highest power of g that divides f is g^v.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x^2*y^2 + 1\n\nflag, q = divides(f*g, f)\nd = divexact(f*g, f)\nv, q = remove(f*g^3, g)\nn = valuation(f*g^3, g)"
},

{
    "location": "mpolynomial_rings/#Ad-hoc-exact-division-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Ad hoc exact division",
    "category": "section",
    "text": "For any ring that implements exact division, the following can be implemented.divexact(f::MyMPoly{T}, c::Integer) where T <: AbstractAlgebra.RingElem\ndivexact(f::MyMPoly{T}, c::Rational) where T <: AbstractAlgebra.RingElem\ndivexact(f::MyMPoly{T}, c::T) where T <: AbstractAlgebra.RingElemDivide the polynomial exactly by the constant c.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = 3x^2*y^2 + 2x + 1\n\nf1 = divexact(f, 5)\nf2 = divexact(f, QQ(2, 3))"
},

{
    "location": "mpolynomial_rings/#Euclidean-division-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Euclidean division",
    "category": "section",
    "text": "Although multivariate polynomial rings are not in general Euclidean, it is possible to define a quotient with remainder function that depends on the polynomial ordering in the case that the quotient ring is a field or a Euclidean domain. In the case that a polynomial g divides a polynomial f, the result no longer depends on the ordering and the remainder is zero, with the quotient agreeing with the exact quotient.divrem(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a tuple (q r) such that f = qg + r, where the coefficients of terms of r whose monomials are divisible by the leading monomial of g are reduced modulo the leading coefficient of g (according to the Euclidean function on the coefficients).Note that the result of this function depends on the ordering of the polynomial ring.div(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemAs per the divrem function, but returning the quotient only. Especially when the quotient happens to be exact, this function can be exceedingly fast.divrem(f::MyMPoly{T}, G::Array{MyMPoly{T}, 1}) where T <: AbstractAlgebra.RingElemAs per the divrem function above, except that each term of r starting with the most significant term, is reduced modulo the leading terms of each of the polynomials in the array G for which the leading monomial is a divisor.A tuple (Q r) is returned from the function, where Q is an array of polynomials of the same length as G, and such that f = r + sum QiGi.The result is again dependent on the ordering in general, but if the polynomials in G are over a field and the reduced generators of a Groebner basis, then the result is unique.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x + y\nh = y + 1\n\nq = div(f, g)\nq, r = divrem(f, g)\nQ, r = divrem(f, [g, h])"
},

{
    "location": "mpolynomial_rings/#Evaluation-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate(f::MyMPoly{T}, A::Array{T, 1}) where T <: AbstractAlgebra.RingElemEvaluate the polynomial f at the values specified by the entries of the array A.evaluate(f::MyMPoly{T}, A::Array{U, 1}) where {T <: AbstractAlgebra.RingElem, U <: Integer}Evaluate the polynomial f at the values specified by the entries of the array A.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\n\nm = evaluate(f, Rational{BigInt}[2, 3])\nn = evaluate(f, [2, 3])"
},

{
    "location": "mpolynomial_rings/#GCD-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "GCD",
    "category": "section",
    "text": "In cases where there is a meaningful Euclidean structure on the coefficient ring, it is possible to compute the GCD of multivariate polynomials.gcd(f::MyMPoly{T}, g::MyMPoly{T}) where T <: AbstractAlgebra.RingElemReturn a greatest common divisor of f and g.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = 2x^2*y + 2x + y + 1\ng = x^2*y^2 + 1\n\nd = gcd(f*g^2, f^2*g)"
},

{
    "location": "mpolynomial_rings/#Interface-for-sparse-distributed,-random-access-multivariates-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Interface for sparse distributed, random access multivariates",
    "category": "section",
    "text": "The following additional functions should be implemented by libraries that provide a sparse distributed polynomial format, stored in a representation for which terms can be accessed in constant time (e.g. where arrays are used to store coefficients and exponent vectors)."
},

{
    "location": "mpolynomial_rings/#Sparse-distributed,-random-access-constructors-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Sparse distributed, random access constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructor, taking arrays of coefficients and exponent vectors, should be provided.(S::MyMPolyRing{T})(A::Vector{T}, m::Vector{Vector{Int}}) where T <: AbstractAlgebra.RingElemCreate the polynomial in the given ring with nonzero coefficients specified by the elements of A and corresponding exponent vectors given by the elements of m.There is no assumption about coefficients being nonzero or terms being in order or unique. Zero terms are removed by the function, duplicate terms are combined (added) and the terms are sorted so that they are in the correct order.Each exponent vector uses a separate integer for each exponent field, the first of which should be the exponent for the most significant variable with respect to the ordering. All exponents must be non-negative.A library may also optionally provide an interface that makes use of BigInt (or any other big integer type) for exponents instead of Int.ExamplesS, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = S(Rational{BigInt}[2, 3, 1], [[3, 2], [1, 0], [0, 1]])"
},

{
    "location": "mpolynomial_rings/#Sparse-distributed,-random-access-basic-manipulation-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Sparse distributed, random access basic manipulation",
    "category": "section",
    "text": "coeff(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the coefficient of the n-th term of f. The first term should be the most significant term with respect to the ordering.coeff(a::MyMPoly{T}, exps::Vector{Int}) where T <: RingElementReturn the coefficient of the term with the given exponent vector, or zero if there is no such term.monomial(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElem\nmonomial!(m::MyMPoly{T}, f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the n-th monomial of f or set m to the n-th monomial of f, respectively. The first monomial should be the most significant term with respect to the ordering. Monomials have coefficient 1 in AbstractAlgebra. See the function term if you also require the coefficient, however, note that only monomials can be compared.term(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemReturn the n-th term of f. The first term should be the one whose monomial is most significant with respect to the ordering.exponent(f::MyMPoly{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElemReturn the exponent of the j-th variable in the n-th term of the polynomial f. The first term is the one with whose monomial is most significant with respect to the ordering.exponent_vector(a::MyMPoly{T}, i::Int) where T <: RingElementReturn a vector of exponents, corresponding to the exponent vector of the i-th term of the polynomial. Term numbering begins at 1 and the exponents are given in the order of the variables for the ring, as supplied when the ring was created.setcoeff!(a::MyMPoly, exps::Vector{Int}, c::S) where S <: RingElementSet the coefficient of the term with the given exponent vector to the given value c. If no such term exists (and c neq 0), one will be inserted. This function takes O(log n) operations if a term with the given exponent already exists and c neq 0, or if the term is inserted at the end of the polynomial. Otherwise it can take O(n) operations in the worst case. This function must return the modified polynomial.ExamplesS, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^3*y + 3x*y^2 + 1\n\nc1 = coeff(f, 1)\nc2 = coeff(f, x^3*y)\nm = monomial(f, 2)\ne1 = exponent(f, 1, 1)\nv1 = exponent_vector(f, 1)\nt1 = term(f, 1)\nsetcoeff!(f, [3, 1], 12)"
},

{
    "location": "mpolynomial_rings/#Unsafe-functions-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Unsafe functions",
    "category": "section",
    "text": "The following functions must be provided, but are considered unsafe, as they may leave the polynomials in an inconsistent state and they mutate their inputs. As usual, such functions should only be applied on polynomials that have no references elsewhere in the system and are mainly intended to be used in carefully written library code, rather than by users.Users should instead build polynomials using the constructors described above.fit!(f::MyMPoly{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial f internally has space for n nonzero terms. This function must mutate the function in-place if it is mutable. It does not return the mutated polynomial. Immutable types can still be supported by defining this function to do nothing.setcoeff!(a::MyMPoly{T}, i::Int, c::T) where T <: RingElement\nsetcoeff!(a::MyMPoly{T}, i::Int, c::U) where {T <: RingElement, U <: Integer}Set the i-th coefficient of the polynomial a to c. No check is performed on the index i or for c = 0. It may be necessary to call combine_like_terms after calls to this function, to remove zero terms. The function must return the modified polynomial.combine_like_terms!(a::MyMPoly{T}) where T <: RingElementRemove zero terms and combine any adjacent terms with the same exponent vector (by adding them). It is assumed that all the exponent vectors are already in the correct order with respect to the ordering. The function must return the resulting polynomial.set_exponent_vector!(a::MyMPoly{T}, i::Int, exps::Vector{Int}) where T <: RingElement Set the i-th exponent vector to the given exponent vector. No check is performed on the index i, which is assumed to be valid (or that the polynomial has enough space allocated). No sorting of exponents is performed by this function. To sort the terms after setting any number of exponents with this function, run the sort_terms! function. The function must return the modified polynomial.sort_terms!(a::MyMPoly{T}) where {T <: RingElement}Sort the terms of the given polynomial according to the polynomial ring ordering. Zero terms and duplicate exponents are ignored. To deal with those call combine_like_terms. The sorted polynomial must be returned by the function."
},

{
    "location": "mpolynomial_rings/#Optional-functionality-for-multivariate-polynomials-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Optional functionality for multivariate polynomials",
    "category": "section",
    "text": "The following functions can optionally be implemented for multivariate polynomial types."
},

{
    "location": "mpolynomial_rings/#Evaluation-2",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Evaluation",
    "category": "section",
    "text": "evaluate(a::MyMPoly{T}, A::Vector{T}) where T <: RingElemEvaluate the polynomial at the given values in the coefficient ring of the polynomial. The result should be an element of the coefficient ring.(a::MyMPoly{T})(vals::Union{NCRingElem, RingElement}...) where T <: RingElementEvaluate the polynomial at the given arguments. This provides functional notation for polynomial evaluation, i.e. f(a b c). It must be defined for each supported polynomial type (Julia does not allow functional notation to be defined for an abstract type).The code for this function in MPoly.jl can be used when implementing this as it provides the most general possible evaluation, which is much more general than the case of evaluation at elements of the same ring.The evaluation should succeed for any set of values for which a multiplication is defined with the product of a coefficient and all the values before it.Note that the values may be in non-commutative rings. Products are performed in the order of the variables in the polynomial ring that the polynomial belongs to, preceded by a multiplication by the coefficient on the left.ExamplesR, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, [\"x\", \"y\"])\nS = MatrixAlgebra(ZZ, 2)\n\nf = x*y + x + y + 1\n\nevaluate(f, [ZZ(1), ZZ(2)])\nf(1, 2)\nf(ZZ(1), ZZ(2))\nf(x - y, x + y)\n\nM1 = S([1 2; 3 4])\nM2 = S([2 4; 1 -1])\nM3 = S([1 -1; 1 1])\n\nf(M1, M2)\nf(M1, ZZ(2))"
},

{
    "location": "mpolynomial_rings/#Derivations-1",
    "page": "Multivariate Polynomial Ring Interface",
    "title": "Derivations",
    "category": "section",
    "text": "The following function allows to compute derivations of multivariate polynomials of type MPoly.derivative(f::MyMPoly{T}, x::MyMPoly{T}) where T <: AbstractAlgebra.RingElemCompute the derivative of f with respect to the variable x.ExamplesR, (x, y) = AbstractAlgebra.PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x*y + x + y + 1\n\nderivative(f, x)\nderivative(f, y)"
},

{
    "location": "mpolynomial/#",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Generic sparse distributed multivariate polynomials",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "mpolynomial/#Generic-sparse-distributed-multivariate-polynomials-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Generic sparse distributed multivariate polynomials",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/MPoly.jl for generic sparse distributed multivariate polynomials over any commutative ring belonging to the AbstractAlgebra abstract type hierarchy.This modules implements the Multivariate Polynomial interface, including the sparse distributed, random access part of the interface.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name.Multivariates are implemented in this module using a Julia array of coefficients and a 2-dimensional Julia array of UInts for the exponent vectors. Note that exponent n is represented by the n-th column of the exponent array, not the n-th row. This is because Julia uses a column major representation."
},

{
    "location": "mpolynomial/#Types-and-parent-objects-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Multivariate polynomials implemented in AbstractAlgebra.jl have type Generic.MPoly{T} where T is the type of elements of the coefficient ring.The polynomials are implemented using a Julia array of coefficients and a 2-dimensional Julia array of UInts for the exponent vectors. Note that exponent n is represented by the n-th column of the exponent array, not the n-th row. This is because Julia uses a column major representation. See the file src/generic/GenericTypes.jl for details.The top bit of each UInt is reserved for overflow detection.Parent objects of such polynomials have type Generic.MPolyRing{T}.The string representation of the variables of the polynomial ring, the base/coefficient ring R and the ordering are stored in the parent object. The polynomial element types belong to the abstract type AbstractAlgebra.MPolyElem{T} and the polynomial ring types belong to the abstract type AbstractAlgebra.MPolyRing{T}.Note that both the generic polynomial ring type Generic.MPolyRing{T} and the abstract type it belongs to, AbstractAlgebra.MPolyRing{T} are both called MPolyRing. The  former is a (parameterised) concrete type for a polynomial ring over a given base ring whose elements have type T. The latter is an abstract type representing all multivariate polynomial ring types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "mpolynomial/#Polynomial-ring-constructors-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Polynomial ring constructors",
    "category": "section",
    "text": "In order to construct multivariate polynomials in AbstractAlgebra.jl, one must first construct the polynomial ring itself. This is accomplished with the following constructor.PolynomialRing(R::AbstractAlgebra.Ring, S::Array{String, 1}; cached::Bool = true, ordering::Symbol=:lex)Given a base ring R and and array S of strings specifying how the generators (variables) should be printed, return a tuple S, (x, ...) representing the new polynomial ring S = Rx ldots and a tuple of the generators (x ) of the ring. By default the parent object S will depend only on R and  (x, ...) and will be cached. Setting the optional argument cached to false will prevent the parent object  S from being cached.The optional named argument ordering can be used to specify an ordering. The currently supported options are :lex, :deglex and :degrevlex.Here are some examples of creating multivariate polynomial rings and making use of the resulting parent objects to coerce various elements into the polynomial ring.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"]; ordering=:deglex)\n\nf = R()\ng = R(123)\nh = R(BigInt(1234))\nk = R(x + 1)\nm = R(x + y + 1)\nderivative(k,x)\nderivative(k,y)All of the examples here are generic polynomial rings, but specialised implementations of polynomial rings provided by external modules will also usually provide a PolynomialRing constructor to allow creation of their polynomial rings."
},

{
    "location": "mpolynomial/#Polynomial-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Polynomial functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": ""
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.vars-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.vars",
    "category": "method",
    "text": "vars(p::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}\n\nReturns the variables actually occuring in p.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.var_index-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.var_index",
    "category": "method",
    "text": "var_index(p::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}\n\nReturn the index of the given variable x. If x is not a variable in a multivariate polynomial ring, an exception is raised.\n\n\n\n"
},

{
    "location": "mpolynomial/#Basic-manipulation-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Basic manipulation",
    "category": "section",
    "text": "vars(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElementvar_index(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElementExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^2 + 2x + 1\n\nV = vars(f)\nvar_index(y) == 2"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.change_base_ring-Union{Tuple{T}, Tuple{MPolyElem{T},Any}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.change_base_ring",
    "category": "method",
    "text": "change_base_ring(p::AbstractAlgebra.MPolyElem{T}, g) where {T <: RingElement}\n\nReturns the polynomial obtained by applying g to the coefficients of p.\n\n\n\n"
},

{
    "location": "mpolynomial/#Changing-base-(coefficient)-rings-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Changing base (coefficient) rings",
    "category": "section",
    "text": "In order to substitute the variables of a polynomial f over a ring T by elements in a T-algebra S, you first have to change the base ring of f using the following function, where g is a function representing the structure homomorphism of the T-algebra S.change_base_ring(p::AbstractAlgebra.MPolyElem{T}, g) where {T <: RingElement}Note that g can also be a Nemo parent, e.g. QQ.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nfz = x^2*y^2 + x + 1\n\nfq = change_base_ring(fz, QQ)"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.coeff-Union{Tuple{T}, Tuple{MPolyElem{T},Array{Int64,1},Array{Int64,1}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.coeff",
    "category": "method",
    "text": "coeff(a::AbstractAlgebra.MPolyElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement\n\nReturn the \"coefficient\" of a (as a multivariate polynomial in the same ring) of the monomial consisting of the product of the variables of the given indices raised to the given exponents (note that not all variables need to appear and the exponents can be zero). E.g. coeff(f, [1, 3], [0, 2]) returns the coefficient of x^0*z^2 in the polynomial f (assuming variables x y z in that order).\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.coeff-Union{Tuple{T}, Tuple{T,Array{T,1},Array{Int64,1}}} where T<:MPolyElem",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.coeff",
    "category": "method",
    "text": "coeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: AbstractAlgebra.MPolyElem\n\nReturn the \"coefficient\" of a (as a multivariate polynomial in the same ring) of the monomial consisting of the product of the given variables to the given exponents (note that not all variables need to appear and the exponents can be zero). E.g. coeff(f, [x, z], [0, 2]) returns the coefficient of x^0*z^2 in the polynomial f. \n\n\n\n"
},

{
    "location": "mpolynomial/#Multivariate-coefficients-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Multivariate coefficients",
    "category": "section",
    "text": "In order to return the \"coefficient\" (as a multivariate polynomial in the same ring), of a given monomial (in which some of the variables may not appear and others may be required to appear to exponent zero), we can use the following function.coeff(a::AbstractAlgebra.MPolyElem{T}, vars::Vector{Int}, exps::Vector{Int}) where T <: RingElement\ncoeff(a::T, vars::Vector{T}, exps::Vector{Int}) where T <: AbstractAlgebra.MPolyElemExamplesR, (x, y, z) = PolynomialRing(ZZ, [\"x\", \"y\", \"z\"])\n\nf = x^4*y^2*z^2 - 2x^4*y*z^2 + 4x^4*z^2 + 2x^2*y^2 + x + 1\n\ncoeff(f, [1, 3], [4, 2]) == coeff(f, [x, z], [4, 2])"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.deflation-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.deflation",
    "category": "method",
    "text": "deflation(f::AbstractAlgebra.MPolyElem{T}) where T <: RingElement\n\nComputes deflation parameters for the exponents of the polynomial f. This is a pair of arrays of integers, the first array of which (the shift) gives the minimum exponent for each variable of the polynomial, and the second of which (the deflation) gives the gcds of all the exponents after subtracting the shift, again per variable. This functionality is used by gcd (and can be used by factorisation algorithms).\n\nend)\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.deflate-Union{Tuple{T}, Tuple{MPolyElem{T},Array{Int64,1},Array{Int64,1}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.deflate",
    "category": "method",
    "text": "deflate(f::AbstractAlgebra.MPolyElem{T}, v::Vector{Int}) where T <: RingElement\n\nReturn a polynomial with the same coefficients as f but whose exponents have been shifted down by the given shifts (supplied as an array of shifts, one for each variable, then deflated (divided) by the given exponents (again supplied as an array of deflation factors, one for each variable). The algorithm automatically replaces a deflation of 0 by 1, to avoid division by 0.  \n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.inflate-Union{Tuple{T}, Tuple{MPolyElem{T},Array{Int64,1},Array{Int64,1}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.inflate",
    "category": "method",
    "text": "inflate(f::AbstractAlgebra.MPolyElem{T}, v::Vector{Int}) where T <: RingElement\n\nReturn a polynomial with the same coefficients as f but whose exponents have been inflated (multiplied) by the given deflation exponents (supplied as an array of inflation factors, one for each variable) and then shifted by the given shifts (again supplied as an array of shifts, one for each variable).  \n\n\n\n"
},

{
    "location": "mpolynomial/#Inflation/deflation-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Inflation/deflation",
    "category": "section",
    "text": "deflation(f::AbstractAlgebra.MPolyElem{T}) where T <: RingElementdeflate(f::AbstractAlgebra.MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElementinflate(f::AbstractAlgebra.MPolyElem{T}, shift::Vector{Int}, defl::Vector{Int}) where T <: RingElementExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = x^7*y^8 + 3*x^4*y^8 - x^4*y^2 + 5x*y^5 - x*y^2\n\ndef, shift = deflation(f)\nf1 = deflate(f, def, shift)\nf2 = inflate(f1, def, shift)\nf2 == f"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.to_univariate-Union{Tuple{T}, Tuple{PolyRing{T},MPolyElem{T}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.to_univariate",
    "category": "method",
    "text": "to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyElem{T}) where T <: AbstractAlgebra.RingElement\n\nAssuming the polynomial p is actually a univariate polynomial, convert the polynomial to a univariate polynomial in the given univariate polynomial ring R. An exception is raised if the polynomial p involves more than one variable.\n\n\n\n"
},

{
    "location": "mpolynomial/#Conversions-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Conversions",
    "category": "section",
    "text": "to_univariate(R::AbstractAlgebra.PolyRing{T}, p::AbstractAlgebra.MPolyElem{T}) where T <: AbstractAlgebra.RingElementExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\nS, z = PolynomialRing(ZZ, \"z\")\n\nf = 2x^5 + 3x^4 - 2x^2 - 1\n\ng = to_univariate(S, f)"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{T}, Tuple{MPolyElem{T},Array{U,1}}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.MPolyElem{T}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}\n\nEvaluate the polynomial by substituting in the array of values for each of the variables. The evaluation will succeed if multiplication is defined between elements of the coefficient ring of a and elements of the supplied vector.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{T}, Tuple{MPolyElem{T},Array{Int64,1},Array{U,1}}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.MPolyElem{T}, vars::Vector{Int}, vals::Vector{U}) where {T <: RingElement, U <: RingElement}\n\nEvaluate the polynomial by substituting in the supplied values in the array vals for the corresponding variables with indices given by the array vars. The evaluation will succeed if multiplication is defined between elements of the coefficient ring of a and elements of vals.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{S}, Tuple{T}, Tuple{S,Array{S,1},Array{U,1}}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where S<:MPolyElem{T} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::S, vars::Vector{S}, vals::Vector{U}) where {S <: AbstractAlgebra.MPolyElem{T}, U <: RingElement} where T <: RingElement\n\nEvaluate the polynomial by substituting in the supplied values in the array vals for the corresponding variables (supplied as polynomials) given by the array vars. The evaluation will succeed if multiplication is defined between elements of the coefficient ring of a and elements of vals.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{T}, Tuple{MPolyElem{T},Array{U,1},Any}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.MPolyElem{T}, A::Vector{U}, g) where {T <: RingElement, U <: RingElement}\n\nEvaluate the polynomial at the supplied values after applying the Map or Function given by g to the coefficients of the polynomial.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{T}, Tuple{MPolyElem{T},Array{Int64,1},Array{U,1},Any}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.MPolyElem{T}, vars::Vector{Int}, vals::Vector{U}, g) where {T <: RingElement, U <: RingElement}\n\nEvaluate the polynomial at the supplied values for the variables with given indices after applying the Map or Function given by g to the coefficients of the polynomial.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{S}, Tuple{T}, Tuple{S,Array{S,1},Array{U,1},Any}} where U<:Union{RingElem, AbstractFloat, Integer, Rational} where S<:MPolyElem{T} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::S, vars::Vector{S}, vals::Vector{U}, g) where {S <: AbstractAlgebra.MPolyElem{T}, U <: RingElement} where T <: RingElement\n\nEvaluate the polynomial at the supplied values for the given variables after  applying the Map or Function given by g to the coefficients of the polynomial.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.evaluate-Union{Tuple{U}, Tuple{T}, Tuple{MPolyElem{T},Array{U,1}}} where U<:NCRingElem where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.evaluate",
    "category": "method",
    "text": "evaluate(a::AbstractAlgebra.MPolyElem{T}, vals::Vector{U}) where {T <: RingElement, U <: NCRingElem}\n\nEvaluate the polynomial at the supplied values, which may be any ring elements, commutative or non-commutative, but in the same ring. Evaluation always proceeds in the order of the variables as supplied when creating the polynomial ring to which a belongs. The evaluation will succeed if a product of a coefficient of the polynomial by one of the values is defined.\n\n\n\n"
},

{
    "location": "mpolynomial/#Evaluation-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Evaluation",
    "category": "section",
    "text": "The following function allows evaluation of a polynomial at all its variables. The result is always in the ring that a product of a coefficient and one of the values belongs to, i.e. if all the values are in the coefficient ring, the result of the evaluation will be too.evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: RingElement}The following functions allow evaluation of a polynomial at some of its variables. Note that the result will be a product of values and an element of the polynomial ring, i.e. even if all the values are in the coefficient ring and all variables are given values, the result will be a constant polynomial, not a coefficient.evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{Int}, ::Vector{U}) where {T <: RingElement, U <: RingElement}evaluate(::S, ::Vector{S}, ::Vector{U}) where {S <: AbstractAlgebra.MPolyElem{T}, U <: RingElement} where T <: RingElementThe following functions allow a map or function to be applied to the coefficients of the polynomial before doing the evaluation. Note that one can simply supply a Ring object for the map (see the examples).evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{U}, ::Any) where {T <: RingElement, U <: RingElement}evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{Int}, ::Vector{U}, ::Any) where {T <: RingElement, U <: RingElement}evaluate(::S, ::Vector{S}, ::Vector{U}, ::Any) where {S <: AbstractAlgebra.MPolyElem{T}, U <: RingElement} where T <: RingElementThe following function allows evaluation of a polynomial at values in a not necessarily commutative ring, e.g. elements of a matrix algebra.evaluate(::AbstractAlgebra.MPolyElem{T}, ::Vector{U}) where {T <: RingElement, U <: NCRingElem}ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nf = 2x^2*y^2 + 3x + y + 1\n\nevaluate(f, BigInt[1, 2])\nevaluate(f, [QQ(1), QQ(2)])\nevaluate(f, [1, 2])\nf(1, 2) == 14\n\nevaluate(f, [x + y, 2y - x])\nf(x + y, 2y - x)\n\nR, (x, y, z) = PolynomialRing(ZZ, [\"x\", \"y\", \"z\"])\n\nf = x^2*y^2 + 2x*z + 3y*z + z + 1\n\nevaluate(f, [1, 3], [3, 4])\nevaluate(f, [x, z], [3, 4])\n\nevaluate(f, [1, 2], [x + z, x - z])\n\nevaluate(f, [2, 4, 6], QQ)\nevaluate(f, [x, z], [2, 4], QQ)\n\nS = MatrixAlgebra(ZZ, 2)\n\nM1 = S([1 2; 3 4])\nM2 = S([2 3; 1 -1])\nM3 = S([-1 1; 1 1])\n\nevaluate(f, [M1, M2, M3])\nf(M1, M2, M3)\n\nf(M1, ZZ(2), M3)\nf(M1, ZZ(2), 3)"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.lc-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.lc",
    "category": "method",
    "text": "lc(p::MPolyElem)\n\nReturn the leading coefficient of the polynomial p.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.lm-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.lm",
    "category": "method",
    "text": "lm(p::MPolyElem)\n\nReturn the leading monomial of the polynomial p.\n\n\n\n"
},

{
    "location": "mpolynomial/#AbstractAlgebra.Generic.lt-Union{Tuple{MPolyElem{T}}, Tuple{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "AbstractAlgebra.Generic.lt",
    "category": "method",
    "text": "lt(p::MPolyElem)\n\nReturn the leading term of the polynomial p.\n\n\n\n"
},

{
    "location": "mpolynomial/#Leading-coefficients,-leading-monomials-and-leading-terms-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Leading coefficients, leading monomials and leading terms",
    "category": "section",
    "text": "The leading coefficient, leading monomial and leading term of a polynomial p are returned by the following functions:lc(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElementlm(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElementlt(p::AbstractAlgebra.MPolyElem{T}) where T <: RingElementExamplesusing AbstractAlgebra\nR,(x,y) = PolynomialRing(ZZ, [\"x\", \"y\"], ordering=:deglex)\np = 2*x*y + 3*y^3\nlt(p)\nlm(p)\nlc(p)\nlt(p) == lc(p) * lm(p)"
},

{
    "location": "mpolynomial/#Base.gcd-Union{Tuple{T}, Tuple{MPoly{T},MPoly{T}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Base.gcd",
    "category": "method",
    "text": "gcd(a::AbstractAlgebra.Generic.MPoly{T}, a::AbstractAlgebra.Generic.MPoly{T}) where {T <: RingElement}\n\nReturn the greatest common divisor of a and b in parent(a).\n\n\n\n"
},

{
    "location": "mpolynomial/#Base.lcm-Union{Tuple{T}, Tuple{MPolyElem{T},MPolyElem{T}}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Base.lcm",
    "category": "method",
    "text": "lcm(a::AbstractAlgebra.MPolyElem{T}, a::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}\n\nReturn the least common multiple of a and b in parent(a).\n\n\n\n"
},

{
    "location": "mpolynomial/#Least-common-multiple,-greatest-common-divisor-1",
    "page": "Generic sparse distributed multivariate polynomials",
    "title": "Least common multiple, greatest common divisor",
    "category": "section",
    "text": "The greated common divisor of two polynomials a and b is returned bygcd(a::AbstractAlgebra.Generic.MPoly{T}, b::AbstractAlgebra.Generic.MPoly{T}) where {T <: RingElement}The least common multiple of two polynomials a and b is returned bylcm(a::AbstractAlgebra.MPolyElem{T}, b::AbstractAlgebra.MPolyElem{T}) where {T <: RingElement}Examplesusing AbstractAlgebra\nR,(x,y) = PolynomialRing(ZZ, [\"x\", \"y\"])\na = x*y + 2*y\nb = x^3*y + y\ngcd(a,b)\nlcm(a,b)\nlcm(a,b) == a * b // gcd(a,b)"
},

{
    "location": "series_rings/#",
    "page": "Series Ring Interface",
    "title": "Series Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "series_rings/#Series-Ring-Interface-1",
    "page": "Series Ring Interface",
    "title": "Series Ring Interface",
    "category": "section",
    "text": "Univariate power series rings are supported in AbstractAlgebra in a variety of different forms, including absolute and relative precision models and Laurent series.In addition to the standard Ring interface, numerous additional functions are required to be present for power series rings."
},

{
    "location": "series_rings/#Types-and-parents-1",
    "page": "Series Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for power series rings and their elements:SeriesRing{T} is the abstract type for all power series ring parent types\nSeriesElem{T} is the abstract type for all power series typesWe have that SeriesRing{T} <: AbstractAlgebra.Ring and  SeriesElem{T} <: AbstractAlgebra.RingElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the coefficient ring of the power series ring. For example, in the case of mathbbZx the type T would be the type of an integer, e.g. BigInt.Within the SeriesElem{T} abstract type is the abstract type RelSeriesElem{T} for relative power series, and AbsSeriesElem{T} for absolute power series.Relative series are typically stored with a valuation and a series that is either zero or that has nonzero constant term. Absolute series are stored starting from the constant term, even if it is zero.If the parent object for a relative series ring over the bignum integers has type MySeriesRing and series in that ring have type MySeries then one would have:MySeriesRing <: SeriesRing{BigInt}\nMySeries <: RelSeriesElem{BigInt}Series rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Series rings should at least be distinguished based on their base (coefficient) ring. But if they have the same base ring and symbol (for their variable/generator) and same default precision, they should certainly have the same parent object.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "series_rings/#Required-functionality-for-series-1",
    "page": "Series Ring Interface",
    "title": "Required functionality for series",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface the Series Ring interface has the following required functions.We suppose that R is a fictitious base ring (coefficient ring) and that S is a series ring over R (e.g. S = Rx) with parent object S of type MySeriesRing{T}. We also assume the series in the ring have type MySeries{T}, where  T is the type of elements of the base (coefficient) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "series_rings/#Constructors-1",
    "page": "Series Ring Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of coefficients, must be available.For relative power series and Laurent series we have:(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int, val::Int) where T <: AbstractAlgebra.RingElemCreate the series in the given ring whose valuation is val, whose absolute precision is given by prec and the coefficients of which are given by A, starting from the first nonzero term. Only len terms of the array are used, the remaining terms being ignored. The value len cannot exceed the length of the supplied array.It is permitted to have trailing zeros in the array, but it is not needed, even if the precision minus the valuation is bigger than the length of the array.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\"; model=:capped_relative)\nT, y = LaurentSeriesRing(ZZ, 10, \"y\")\nU, z = LaurentSeriesField(QQ, 10, \"z\")\n \nf = S(Rational{BigInt}[2, 3, 1], 3, 6, 2)\ng = T(BigInt[2, 3, 1], 3, 6, 2)\nh = U(Rational{BigInt}[2, 3, 1], 3, 6, 2)For absolute power series we have:(S::MySeriesRing{T})(A::Array{T, 1}, len::Int, prec::Int) where T <: AbstractAlgebra.RingElemCreate the series in the given ring whose absolute precision is given by prec and the coefficients of which are given by A, starting from the constant term. Only len terms of the array are used, the remaining terms being ignored.Note that len is usually maintained separately of any polynomial that is underlying the power series. This allows for easy trucation of a power series without actually modifying the polynomial underlying it.It is permitted to have trailing zeros in the array, but it is not needed, even if the precision is bigger than the length of the array.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\"; model=:capped_absolute)\n\nf = S(Rational{BigInt}[0, 2, 3, 1], 4, 6)"
},

{
    "location": "series_rings/#Data-type-and-parent-object-methods-1",
    "page": "Series Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "var(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn a Symbol representing the variable (generator) of the series ring. Note that this is a Symbol not a String, though its string value will usually be used when printing series.max_precision(S::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn the (default) maximum precision of the power series ring. This is the precision that the output of an operation will be if it cannot be represented to full precision (e.g. because it mathematically has infinite precision).This value is usually supplied upon creation of the series ring and stored in the ring. It is independent of the precision which each series in the ring actually has. Those are stored on a per element basis in the actual series elements.ExamplesS, x = PowerSeriesRing(QQ, 10, \"x\")\n\nvsym = var(S)\nmax_precision(S) == 10"
},

{
    "location": "series_rings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Series Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "pol_length(f::MySeries{T}) where T <: AbstractAlgebra.RingElemReturn the length of the polynomial underlying the given power series. This is not generally useful to the user, but is used internally.set_length!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElemThis function sets the effective length of the polynomial underlying the given series. The function doesn\'t modify the actual polynomial, but simply changes the number of terms of the polynomial which are considered to belong to the power series. The remaining terms are ignored.This function cannot set the length to a value greater than the length of any underlying polynomial.The function mutates the series in-place but does not return the mutated series.precision(f::MySeries{T})Returns the absolute precision of f.set_prec!(f::MySeries{T}, prec::Int)Set the absolute precision of the given series to the given value.This function mutates the series in-place but does not return the mutated series.valuation(f::MySeries{T})Return the valuation of the given series.set_val!(f::MySeries{T}, val::Int)For relative series and Laurent series only, this function alters the valuation of the given series to the given value.The series is mutated in-place but does not return the mutated series.polcoeff(f::MySeries{T}, n::Int) Return the coefficient of degree n of the polynomial underlying the series. If n is larger than the degree of this polynomial, zero is returned. This function is not generally of use to the user but is used internally.setcoeff!(f::MySeries{T}, n::Int, a::T) where T <: AbstractAlgebra.RingElemSet the degree n coefficient of the polynomial underlying f to a. This mutates the polynomial in-place if possible and returns the mutated series (so that immutable types can also be supported). The function must not assume that the polynomial already has space for n + 1 coefficients. The polynomial must be resized if this is not the case.Note that this function is not required to normalise the polynomial and is not necessarily useful to the user, but is used extensively by the generic functionality in AbstractAlgebra.jl. It is for setting raw coefficients in the representation.normalise(f::MySeries{T}, n::Int)Given a series f represented by a polynomial of at least the given length, return the normalised length of the underlying polynomial assuming it has length at most n. This function does not actually normalise the polynomial and is not particularly useful to the user. It is used internally.renormalize!(f::MySeries{T}) where T <: AbstractAlgebra.RingElemGiven a relative series or Laurent series whose underlying polynomial has zero constant term, say as the result of some internal computation, renormalise the series so that the  polynomial has nonzero constant term. The precision and valuation of the series are adjusted to compensate. This function is not intended to be useful to the user, but is  used internally.fit!(f::MySeries{T}, n::Int) where T <: AbstractAlgebra.RingElemEnsure that the polynomial underlying f internally has space for n coefficients. This function must mutate the series in-place if it is mutable. It does not return the mutated series. Immutable types can still be supported by defining this function to do nothing.Some interfaces for C polynomial types automatically manage the internal allocation of polynomials in every function that can be called on them. Explicit adjustment by the generic code in AbstractAlgebra.jl is not required. In such cases, this function can also be defined to do nothing.gen(R::MySeriesRing{T}) where T <: AbstractAlgebra.RingElemReturn the generator x of the series ring.ExamplesS, x = PowerSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 3x + x^3 + O(x^5)\ng = S(BigInt[1, 2, 0, 1, 0, 0, 0], 4, 10, 3);\n\nn = pol_length(f)\nc = polcoeff(f, 1)\nset_length!(g, 3)\ng = setcoeff!(g, 2, BigInt(11))\nfit!(g, 8)\ng = setcoeff!(g, 7, BigInt(4))\nw = gen(S)\nisgen(w) == true"
},

{
    "location": "series/#",
    "page": "Generic power series",
    "title": "Generic power series",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "series/#Generic-power-series-1",
    "page": "Generic power series",
    "title": "Generic power series",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of capped relative and absolute power series over  any computable commutative ring R.Capped relative power series are power series of the form a_jx^j + a_j+1x^j+1 + cdots + a_k-1x^k-1 + O(x^k) where a_j in R and the relative precision k - j is at most equal to some specified precision n.Capped absolute power series are power series of the form a_jx^j + a_j+1x^j+1 + cdots + a_n-1x^n-1 + O(x^n) where j geq 0, a_j in R and the precision n is fixed.There are two implementations of relative series: relative power series, implemented in src/generic/RelSeries.jl for which j  0 in the above description, and Laurent series where j can be negative, implemented in src/generic/Laurent.jl. Note that there are two implementations for Laurent series, one over rings and one over fields, though in practice most of the implementation uses the same code in both cases.There is a single implementation of absolute series: absolute power series, implemented in src/generic/AbsSeries.jl.As well as implementing the Series Ring interface, the series modules in AbstractAlgebra.jl implement the generic algorithms described below.All of the generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default so that it is not necessary to qualify function names."
},

{
    "location": "series/#Types-and-parent-objects-1",
    "page": "Generic power series",
    "title": "Types and parent objects",
    "category": "section",
    "text": "The types of generic series implemented by AbstractAlgebra.jl are Generic.RelSeries{T}, Generic.AbsSeries{T}, Generic.LaurentSeriesRingElem{T} and Generic.LaurentSeriesFieldElem{T}.Relative power series elements belong to the abstract type AbstractAlgebra.RelSeriesElem.Laurent series elements belong directly to either AbstractAlgebra.RingElem or AbstractAlgebra.FieldElem since it is more useful to be able to distinguish whether they belong to a ring or field than it is to distinguish that they are relative series.Absolute power series elements belong to AbstractAlgebra.AbsSeriesElem.The parent types for relative and absolute power series, Generic.RelSeriesRing{T}  and Generic.AbsSeriesRing{T} respectively, belong to AbstractAlgebra.SeriesRing{T}.The parent types for Laurent series rings and fields, Generic.LaurentSeriesRing{T} and Generic.LaurentSeriesField{T} respectively, belong directly to  AbstractAlgebra.Ring and AbstractAlgebra.Field respectively.The default precision, string representation of the variable and base ring R of a generic power series are stored in its parent object. "
},

{
    "location": "series/#Series-ring-constructors-1",
    "page": "Generic power series",
    "title": "Series ring constructors",
    "category": "section",
    "text": "In order to construct series in AbstractAlgebra.jl, one must first construct the ring itself. This is accomplished with any of the following constructors.PowerSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true, model=:capped_relative)LaurentSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true)LaurentSeriesRing(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)Given a base ring R, a maximum precision (relative or absolute, depending on the model) and a string s specifying how the generator (variable) should be printed, return a typle S, x representing the series ring and its generator.By default, S will depend only on S, x and the maximum precision and will be cached. Setting the optional argument cached to false will prevent this.In the case of power series, the optional argument model can be set to either :capped_absolute or capped_relative, depending on which power series model is required.Here are some examples of constructing various kinds of series rings and coercing various elements into those rings.ExamplesR, x = PowerSeriesRing(ZZ, 10, \"x\")\nS, y = PowerSeriesRing(ZZ, 10, \"y\"; model=:capped_absolute)\nT, z = LaurentSeriesRing(ZZ, 10, \"z\")\nU, w = LaurentSeriesField(QQ, 10, \"w\")\n\nf = R()\ng = S(123)\nh = U(BigInt(1234))\nk = T(z + 1)"
},

{
    "location": "series/#Big-oh-notation-1",
    "page": "Generic power series",
    "title": "Big-oh notation",
    "category": "section",
    "text": "Series elements can be given a precision using the big-oh notation. This is provided by a function of the following form, (or something equivalent for Laurent series):O(x::SeriesElem)ExamplesR, x = PowerSeriesRing(ZZ, 10, \"x\")\nS, y = LaurentSeriesRing(ZZ, 10, \"y\")\n\nf = 1 + 2x + O(x^5)\ng = 2y + 7y^2 + O(y^7)What is happening here in practice is that O(x^n) is creating the series 0 + O(x^n) and the rules for addition of series dictate that if this is added to a series of  greater precision, then the lower of the two precisions must be used.Of course it may be that the precision of the series that O(x^n) is added to is already lower than n, in which case adding O(x^n) has no effect. This is the case if the default precision is too low, since x on its own has the default precision."
},

{
    "location": "series/#Power-series-models-1",
    "page": "Generic power series",
    "title": "Power series models",
    "category": "section",
    "text": "Capped relative power series have their maximum relative precision capped at some value prec_max. This means that if the leading term of a nonzero power series element is c_ax^a and the precision is b then the power series is of the form  c_ax^a + c_a+1x^a+1 + ldots + O(x^a + b).The zero power series is simply taken to be 0 + O(x^b).The capped relative model has the advantage that power series are stable multiplicatively. In other words, for nonzero power series f and g we have that divexact(f*g), g) == f.However, capped relative power series are not additively stable, i.e. we do not always have (f + g) - g = f.Similar comments apply to Laurent series.On the other hand, capped absolute power series have their absolute precision capped. This means that if the leading term of a nonzero power series element is c_ax^a and the precision is b then the power series is of the form c_ax^a + c_a+1x^a+1 + ldots + O(x^b).Capped absolute series are additively stable, but not necessarily multiplicatively stable.For all models, the maximum precision is also used as a default precision in the case of coercing coefficients into the ring and for any computation where the result could mathematically be given to infinite precision.In all models we say that two power series are equal if they agree up to the minimum absolute precision of the two power series.Thus, for example, x^5 + O(x^10) == 0 + O(x^5), since the minimum absolute precision is 5.During computations, it is possible for power series to lose relative precision due to cancellation. For example if f = x^3 + x^5 + O(x^8) and g = x^3 + x^6 + O(x^8) then f - g = x^5 - x^6 + O(x^8) which now has relative precision 3 instead of relative precision 5.Amongst other things, this means that equality is not transitive. For example x^6 + O(x^11) == 0 + O(x^5) and x^7 + O(x^12) == 0 + O(x^5) but x^6 + O(x^11) neq x^7 + O(x^12).Sometimes it is necessary to compare power series not just for arithmetic equality, as above, but to see if they have precisely the same precision and terms. For this purpose we introduce the isequal function.For example, if f = x^2 + O(x^7) and g = x^2 + O(x^8) and h = 0 + O(x^2) then f == g, f == h and g == h, but isequal(f, g), isequal(f, h) and isequal(g, h) would all return false. However, if k = x^2 + O(x^7) then isequal(f, k) would return true.There are further difficulties if we construct polynomial over power series. For example, consider the polynomial in y over the power series ring in x over the rationals. Normalisation of such polynomials is problematic. For instance, what is the leading coefficient of (0 + O(x^10))y + (1 + O(x^10))?If one takes it to be (0 + O(x^10)) then some functions may not terminate due to the fact that algorithms may require the degree of polynomials to decrease with each iteration. Instead, the degree may remain constant and simply accumulate leading terms which are arithmetically zero but not identically zero.On the other hand, when constructing power series over other power series, if we simply throw away terms which are arithmetically equal to zero, our computations may have different output depending on the order in which the power series are added!One should be aware of these difficulties when working with power series. Power series, as represented on a computer, simply don\'t satisfy the axioms of a ring. They must be used with care in order to approximate operations in a mathematical power series ring.Simply increasing the precision will not necessarily give a \"more correct\" answer and some computations may not even terminate due to the presence of arithmetic zeroes!An absolute power series ring over a ring R with precision p behaves  very much like the quotient Rx(x^p) of the polynomial ring over R. Therefore one can often treat absolute power series rings as though they were rings. However, this depends on all series being given a precision equal to the specified maximum precision and not a lower precision."
},

{
    "location": "series/#Basic-ring-functionality-1",
    "page": "Generic power series",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "All power series models provide the functionality described in the Ring and Series Ring interfaces.ExamplesS, x = PowerSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 3x + x^3 + O(x^10)\ng = 1 + 2x + x^2 + O(x^10)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nn = pol_length(f)\nc = polcoeff(f, 3)\nU = base_ring(S)\nv = var(S)\nT = parent(x + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)\np = precision(f)"
},

{
    "location": "series/#Series-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic power series",
    "title": "Series functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below is automatically provided by AbstractAlgebra.jl for any series module that implements the full Series Ring interface. This includes AbstractAlgebra\'s own generic series rings.Of course, modules are encouraged to provide specific implementations of the functions described here, that override the generic implementation.Unless otherwise noted, the functions are available for all series models, including Laurent series. We denote this by using the abstract type AbstractAlgebra.RelSeriesElem, even though absolute series and Laurent series types do not belong to this abstract type."
},

{
    "location": "series/#AbstractAlgebra.Generic.modulus-Union{Tuple{SeriesElem{T}}, Tuple{T}} where T<:ResElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "method",
    "text": "modulus(a::AbstractAlgebra.SeriesElem{T}) where {T <: ResElem}\n\nReturn the modulus of the coefficients of the given power series.\n\n\n\n"
},

{
    "location": "series/#AbstractAlgebra.Generic.isgen-Tuple{RelSeriesElem}",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "method",
    "text": "isgen(a::AbstractAlgebra.RelSeriesElem)\n\nReturn true if the given power series is arithmetically equal to the generator of its power series ring to its current precision, otherwise return false.\n\n\n\n"
},

{
    "location": "series/#AbstractAlgebra.Generic.isunit-Tuple{RelSeriesElem}",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::AbstractAlgebra.RelSeriesElem)\n\nReturn true if the given power series is arithmetically equal to a unit, i.e. is invertible, otherwise return false.\n\n\n\n"
},

{
    "location": "series/#Basic-functionality-1",
    "page": "Generic power series",
    "title": "Basic functionality",
    "category": "section",
    "text": "coeff(a::AbstractAlgebra.SeriesElem, n::Int)Return the degree n coefficient of the given power series. Note coefficients are numbered from n = 0 for the constant coefficient. If n exceeds the current precision of the power series, the function returns a zero coefficient.For power series types, n must be non-negative. Laurent series do not have this restriction.modulus{T <: ResElem}(::SeriesElem{T})isgen(::RelSeriesElem)isunit(::RelSeriesElem)ExamplesR, t = PowerSeriesRing(QQ, 10, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = O(x^4)\nb = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)\n\nk = isgen(gen(R))\nm = isunit(-1 + x + 2x^2)\nn = valuation(a)\np = valuation(b)\nc = coeff(b, 2)"
},

{
    "location": "series/#AbstractAlgebra.Generic.shift_left-Union{Tuple{T}, Tuple{RelSeriesElem{T},Int64}} where T<:RingElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.shift_left",
    "category": "method",
    "text": "shift_left(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}\n\nReturn the power series x shifted left by n terms, i.e. multiplied by x^n.\n\n\n\n"
},

{
    "location": "series/#AbstractAlgebra.Generic.shift_right-Union{Tuple{T}, Tuple{RelSeriesElem{T},Int64}} where T<:RingElem",
    "page": "Generic power series",
    "title": "AbstractAlgebra.Generic.shift_right",
    "category": "method",
    "text": "shift_right(x::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}\n\nReturn the power series x shifted right by n terms, i.e. divided by x^n.\n\n\n\n"
},

{
    "location": "series/#Shifting-1",
    "page": "Generic power series",
    "title": "Shifting",
    "category": "section",
    "text": "shift_left{T <: RingElem}(::RelSeriesElem{T}, ::Int)shift_right{T <: RingElem}(::RelSeriesElem{T}, ::Int)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 2x + x^3\nb = O(x^4)\nc = 1 + x + 2x^2 + O(x^5)\nd = 2x + x^3 + O(x^4)\n\nf = shift_left(a, 2)\ng = shift_left(b, 2)\nh = shift_right(c, 1)\nk = shift_right(d, 3)"
},

{
    "location": "series/#Base.truncate-Union{Tuple{T}, Tuple{RelSeriesElem{T},Int64}} where T<:RingElem",
    "page": "Generic power series",
    "title": "Base.truncate",
    "category": "method",
    "text": "truncate(a::AbstractAlgebra.RelSeriesElem{T}, n::Int) where {T <: RingElement}\n\nReturn a truncated to (absolute) precision n.\n\n\n\n"
},

{
    "location": "series/#Truncation-1",
    "page": "Generic power series",
    "title": "Truncation",
    "category": "section",
    "text": "truncate{T <: RingElem}(::RelSeriesElem{T}, ::Int)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 2x + x^3\nb = O(x^4)\nc = 1 + x + 2x^2 + O(x^5)\nd = 2x + x^3 + O(x^4)\n\nf = truncate(a, 3)\ng = truncate(b, 2)\nh = truncate(c, 7)\nk = truncate(d, 5)"
},

{
    "location": "series/#Base.inv-Tuple{RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the inverse of the power series a, i.e. 1a.\n\n\n\n"
},

{
    "location": "series/#Division-1",
    "page": "Generic power series",
    "title": "Division",
    "category": "section",
    "text": "inv(::RelSeriesElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\n\na = 1 + x + 2x^2 + O(x^5)\nb = S(-1)\n\nc = inv(a)\nd = inv(b)"
},

{
    "location": "series/#Base.exp-Tuple{RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.exp",
    "category": "method",
    "text": "exp(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the exponential of the power series a.\n\n\n\n"
},

{
    "location": "series/#Base.sqrt-Tuple{RelSeriesElem}",
    "page": "Generic power series",
    "title": "Base.sqrt",
    "category": "method",
    "text": "sqrt(a::AbstractAlgebra.RelSeriesElem)\n\nReturn the square root of the power series a.\n\n\n\n"
},

{
    "location": "series/#Special-functions-1",
    "page": "Generic power series",
    "title": "Special functions",
    "category": "section",
    "text": "Base.exp(a::RelSeriesElem)Base.sqrt(a::RelSeriesElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PowerSeriesRing(R, 30, \"x\")\nT, z = PowerSeriesRing(QQ, 30, \"z\")\n\na = 1 + z + 3z^2 + O(z^5)\nb = z + 2z^2 + 5z^3 + O(z^5)\n\nc = exp(x + O(x^40))\nd = divexact(x, exp(x + O(x^40)) - 1)\nf = exp(b)\nh = sqrt(a)"
},

{
    "location": "puiseux/#",
    "page": "Generic Puiseux series",
    "title": "Generic Puiseux series",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "puiseux/#Generic-Puiseux-series-1",
    "page": "Generic Puiseux series",
    "title": "Generic Puiseux series",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of Puiseux series over any computable commutative ring R.Puiseux series are power series of the form a_jx^jm + a_j+1x^(j+1)m + cdots + a_k-1x^(k-1)m + O(x^km) for some integer m  0 where i geq 0, a_i in R and the relative precision k - j is at most equal to some specified precision n.The generic Puiseux series module is implemented in src/generic/PuiseuxSeries.jl.As well as implementing the Series Ring interface, the Puiseux series module in AbstractAlgebra.jl implements the generic algorithms described below.All of the generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default so that it is not necessary to qualify function names."
},

{
    "location": "puiseux/#Types-and-parent-objects-1",
    "page": "Generic Puiseux series",
    "title": "Types and parent objects",
    "category": "section",
    "text": "The types of generic polynomials implemented by AbstractAlgebra.jl are Generic.PuiseuxSeriesRingElem{T} and Generic.PuiseuxSeriesFieldElem{T}.Both series element types belong to the union type Generic.PuiseuxSeriesElem.Puiseux series elements belong directly to either AbstractAlgebra.RingElem or AbstractAlgebra.FieldElem since it is more useful to be able to distinguish whether they belong to a ring or field than it is to distinguish that they are Puiseux series.The parent types for Puiseux series, Generic.PuiseuxSeriesRing{T}  and Generic.PuiseuxSeriesField{T} respectively, belong to AbstractAlgebra.Ring and AbstractAlgebra.Field respectively.The default precision, string representation of the variable and base ring R of a generic Puiseux series are stored in its parent object. "
},

{
    "location": "puiseux/#Puisuex-series-ring-constructors-1",
    "page": "Generic Puiseux series",
    "title": "Puisuex series ring constructors",
    "category": "section",
    "text": "In order to construct Puiseux series in AbstractAlgebra.jl, one must first construct the ring itself. This is accomplished with any of the following constructors.PuiseuxSeriesRing(R::AbstractAlgebra.Ring, prec_max::Int, s::AbstractString; cached::Bool = true)PuiseuxSeriesRing(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)PuiseuxSeriesField(R::AbstractAlgebra.Field, prec_max::Int, s::AbstractString; cached::Bool = true)Given a base ring R, a maximum relative precision and a string s specifying how the generator (variable) should be printed, return a typle S, x representing the Puiseux series ring and its generator.By default, S will depend only on S, x and the maximum precision and will be cached. Setting the optional argument cached to false will prevent this.Here are some examples of constructing various kinds of Puiseux series rings and coercing various elements into those rings.ExamplesR, x = PuiseuxSeriesRing(ZZ, 10, \"x\")\nS, y = PuiseuxSeriesField(QQ, 10, \"y\")\n\nf = R()\ng = S(123)\nh = R(BigInt(1234))\nk = S(y + 1)"
},

{
    "location": "puiseux/#Big-oh-notation-1",
    "page": "Generic Puiseux series",
    "title": "Big-oh notation",
    "category": "section",
    "text": "Series elements can be given a precision using the big-oh notation. This is provided by a function of the following form, (or something equivalent for Laurent series):O(x::SeriesElem)ExamplesR, x = PuiseuxSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 2x + O(x^5)\ng = 2x^(1//3) + 7y^(2//3) + O(y^(7//3))What is happening here in practice is that O(x^n) is creating the series 0 + O(x^n) and the rules for addition of series dictate that if this is added to a series of  greater precision, then the lower of the two precisions must be used.Of course it may be that the precision of the series that O(x^n) is added to is already lower than n, in which case adding O(x^n) has no effect. This is the case if the default precision is too low, since x on its own has the default precision."
},

{
    "location": "puiseux/#Puiseux-series-implementation-1",
    "page": "Generic Puiseux series",
    "title": "Puiseux series implementation",
    "category": "section",
    "text": "Puiseux series have their maximum relative precision capped at some value prec_max. This refers to the internal Laurent series used to store the Puiseux series, i.e. the series without denominators in the exponents.The Puiseux series type stores such a Laurent series and a scale or denominator for the exponents. For example, f(x) = 1 + x^13 + 2x^23 + O(x^73) would be stored as a Laurent series 1 + x + 2x^2 + O(x^7) and a scale of 3..The maximum precision is also used as a default (Laurent) precision in the case of coercing coefficients into the ring and for any computation where the result could mathematically be given to infinite precision.In all models we say that two Puiseux series are equal if they agree up to the minimum absolute precision of the two power series.Thus, for example, x^5 + O(x^10) == 0 + O(x^5), since the minimum absolute precision is 5.Sometimes it is necessary to compare Puiseux series not just for arithmetic equality, as above, but to see if they have precisely the same precision and terms. For this purpose we introduce the isequal function.For example, if f = x^2 + O(x^7) and g = x^2 + O(x^8) and h = 0 + O(x^2) then f == g, f == h and g == h, but isequal(f, g), isequal(f, h) and isequal(g, h) would all return false. However, if k = x^2 + O(x^7) then isequal(f, k) would return true.There are a number of technicalities that must be observed when working with Puiseux series. As these are the same as for the other series rings in AbstractAlgebra.jl, we refer the reader to the documentation of series rings for information about these issues."
},

{
    "location": "puiseux/#Basic-ring-functionality-1",
    "page": "Generic Puiseux series",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "All Puiseux series provide the functionality described in the Ring and Series Ring interfaces with the exception of the pol_length and polcoeff functions.ExamplesS, x = PuiseuxSeriesRing(ZZ, 10, \"x\")\n\nf = 1 + 3x + x^3 + O(x^10)\ng = 1 + 2x^(1//3) + x^(2//3) + O(x^(7//3))\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nU = base_ring(S)\nv = var(S)\nT = parent(x + 1)\ng == deepcopy(g)\nt = divexact(2g, 2)\np = precision(f)"
},

{
    "location": "puiseux/#Puiseux-series-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic Puiseux series",
    "title": "Puiseux series functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below is automatically provided by AbstractAlgebra.jl for any Puiseux series.Of course, modules are encouraged to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "puiseux/#AbstractAlgebra.Generic.modulus-Union{Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}}}, Tuple{T}} where T<:ResElem",
    "page": "Generic Puiseux series",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "method",
    "text": "modulus(a::Generic.PuiseuxSeriesElem{T}) where {T <: ResElem}\n\nReturn the modulus of the coefficients of the given Puiseux series.\n\n\n\n"
},

{
    "location": "puiseux/#AbstractAlgebra.Generic.isgen-Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}}",
    "page": "Generic Puiseux series",
    "title": "AbstractAlgebra.Generic.isgen",
    "category": "method",
    "text": "isgen(a::Generic.PuiseuxSeriesElem)\n\nReturn true if the given Puiseux series is arithmetically equal to the generator of its Puiseux series ring to its current precision, otherwise return false.\n\n\n\n"
},

{
    "location": "puiseux/#AbstractAlgebra.Generic.isunit-Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}}",
    "page": "Generic Puiseux series",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::Generic.PuiseuxSeriesElem)\n\nReturn true if the given Puiseux series is arithmetically equal to a unit, i.e. is invertible, otherwise return false.\n\n\n\n"
},

{
    "location": "puiseux/#Basic-functionality-1",
    "page": "Generic Puiseux series",
    "title": "Basic functionality",
    "category": "section",
    "text": "coeff(a::Generic.PuiseuxSeriesElem, n::Int)coeff(a::Generic.PuiseuxSeriesElem, n::Rational{Int})Return the coefficient of the term of exponent n of the given power series. If n exceeds the current precision of the power series or does not correspond to a nonzero term of the Puiseux series, the function returns a zero coefficient.modulus{T <: ResElem}(::Generic.PuiseuxSeriesElem{T})isgen(::Generic.PuiseuxSeriesElem)isunit(::Generic.PuiseuxSeriesElem)ExamplesR, t = PuiseuxSeriesRing(QQ, 10, \"t\")\nS, x = PuiseuxSeriesRing(R, 30, \"x\")\n\na = O(x^4)\nb = (t + 3)*x + (t^2 + 1)*x^2 + O(x^4)\n\nk = isgen(gen(R))\nm = isunit(-1 + x^(1//3) + 2x^2)\nn = valuation(a)\np = valuation(b)\nc = coeff(b, 2)"
},

{
    "location": "puiseux/#Base.inv-Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}}",
    "page": "Generic Puiseux series",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nGiven a non-singular ntimes n matrix over a ring the tuple X d consisting of an ntimes n matrix X and a denominator d such that AX = dI_n, where I_n is the ntimes n identity matrix. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "puiseux/#Division-1",
    "page": "Generic Puiseux series",
    "title": "Division",
    "category": "section",
    "text": "inv(::Generic.PuiseuxSeriesElem)ExamplesR, x = PuiseuxSeriesRing(QQ, 30, \"x\")\n\na = 1 + x + 2x^2 + O(x^5)\nb = R(-1)\n\nc = inv(a)\nd = inv(b)"
},

{
    "location": "puiseux/#Base.exp-Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}}",
    "page": "Generic Puiseux series",
    "title": "Base.exp",
    "category": "method",
    "text": "exp(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement\n\nReturn the exponential of the given Puiseux series a.\n\n\n\n"
},

{
    "location": "puiseux/#Base.sqrt-Tuple{Union{PuiseuxSeriesFieldElem{T}, PuiseuxSeriesRingElem{T}} where T<:Union{RingElem, AbstractFloat, Integer, Rational}}",
    "page": "Generic Puiseux series",
    "title": "Base.sqrt",
    "category": "method",
    "text": "sqrt(a::Generic.PuiseuxSeriesElem{T}) where T <: RingElement\n\nReturn the square root of the given Puiseux series a.\n\n\n\n"
},

{
    "location": "puiseux/#Special-functions-1",
    "page": "Generic Puiseux series",
    "title": "Special functions",
    "category": "section",
    "text": "Base.exp(a::Generic.PuiseuxSeriesElem)Base.sqrt(a::Generic.PuiseuxSeriesElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS, x = PuiseuxSeriesRing(R, 30, \"x\")\nT, z = PuiseuxSeriesRing(QQ, 30, \"z\")\n\na = 1 + z + 3z^2 + O(z^5)\nb = z + 2z^2 + 5z^3 + O(z^5)\n\nc = exp(x + O(x^40))\nd = divexact(x, exp(x + O(x^40)) - 1)\nf = exp(b)\nh = sqrt(a)"
},

{
    "location": "residue_rings/#",
    "page": "Residue Ring Interface",
    "title": "Residue Ring Interface",
    "category": "page",
    "text": ""
},

{
    "location": "residue_rings/#Residue-Ring-Interface-1",
    "page": "Residue Ring Interface",
    "title": "Residue Ring Interface",
    "category": "section",
    "text": "Residue rings (currently a quotient ring modulo a principal ideal) are supported in AbstractAlgebra.jl, at least for Euclidean base rings. In addition to the standard Ring interface, some additional functions are required to be present for residue rings."
},

{
    "location": "residue_rings/#Types-and-parents-1",
    "page": "Residue Ring Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides four abstract types for residue rings and their elements:ResRing{T} is the abstract type for residue ring parent types\nResField{T} is the abstract type for residue rings known to be fields\nResElem{T} is the abstract type for types of elements of residue rings (residues)\nResFieldElem{T} is the abstract type for types of elements of residue fieldsWe have that ResRing{T} <: AbstractAlgebra.Ring and  ResElem{T} <: AbstractAlgebra.RingElem.Note that these abstract types are parameterised. The type T should usually be the type of elements of the base ring of the residue ring/field.If the parent object for a residue ring has type MyResRing and residues in that ring have type MyRes then one would have:MyResRing <: ResRing{BigInt}\nMyRes <: ResElem{BigInt}Residue rings should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Residue rings should at least be distinguished based on their base ring and modulus (the principal ideal one is taking a quotient of the base ring by).See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "residue_rings/#Required-functionality-for-residue-rings-1",
    "page": "Residue Ring Interface",
    "title": "Required functionality for residue rings",
    "category": "section",
    "text": "In addition to the required functionality for the Ring interface the Residue Ring interface has the following required functions.We suppose that R is a fictitious base ring, m is an element of that ring, and that S is the residue ring (quotient ring) R(m) with parent object S of type MyResRing{T}. We also assume the residues r pmodm in the residue ring have type MyRes{T}, where T is the type of elements of the base ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "residue_rings/#Data-type-and-parent-object-methods-1",
    "page": "Residue Ring Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "modulus(S::MyResRing{T}) where T <: AbstractAlgebra.RingElemReturn the modulus of the given residue ring, i.e. if the residue ring S was specified to be R(m), return m.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nm = modulus(S)"
},

{
    "location": "residue_rings/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Residue Ring Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "data(f::MyRes{T}) where T <: AbstractAlgebra.RingElemGiven a residue r pmodm, represented as such, return r.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S(x^2 + 2)\n\nd = data(f)"
},

{
    "location": "residue/#",
    "page": "Generic residue rings",
    "title": "Generic residue rings",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "residue/#Generic-residue-rings-1",
    "page": "Generic residue rings",
    "title": "Generic residue rings",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Residue.jl for generic residue rings over any Euclidean domain (in practice most of the functionality is provided for GCD domains that provide a meaningful GCD function) belonging to the AbstractAlgebra.jl abstract type hierarchy.As well as implementing the Residue Ring interface a number of generic algorithms are implemented for residue rings. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "residue/#Types-and-parent-objects-1",
    "page": "Generic residue rings",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Residues implemented using the AbstractAlgebra generics have type Generic.Res{T} or in the case of residue rings that are known to be fields, Generic.ResF{T}, where T is the type of elements of the base ring. See the file src/generic/GenericTypes.jl for details.Parent objects of residue ring elements have type Generic.ResRing{T} and those of residue fields have type GenericResField{T}.The defining modulus of the residue ring is stored in the parent object.The residue element types belong to the abstract type AbstractAlgebra.ResElem{T} or AbstractAlgebra.ResFieldElem{T} in the case of residue fields, and the residue ring types belong to the abstract type AbstractAlgebra.ResRing{T} or AbstractAlgebra.ResField{T} respectively. This enables one to write generic functions that can accept any AbstractAlgebra residue type.Note that both the generic residue ring type Generic.ResRing{T} and the abstract type it belongs to, AbstractAlgebra.ResRing{T} are both called ResRing, and  similarly for the residue field types. In each case, the  former is a (parameterised) concrete type for a residue ring over a given base ring whose elements have type T. The latter is an abstract type representing all residue ring types in  AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "residue/#Residue-ring-constructors-1",
    "page": "Generic residue rings",
    "title": "Residue ring constructors",
    "category": "section",
    "text": "In order to construct residues in AbstractAlgebra.jl, one must first construct the resiude ring itself. This is accomplished with one of the following constructors.ResidueRing(R::AbstractAlgebra.Ring, m::AbstractAlgebra.RingElem; cached::Bool = true)ResidueField(R::AbstractAlgebra.Ring, m::AbstractAlgebra.RingElem; cached::Bool = true)Given a base ring R and residue m contained in this ring, return the parent object of the residue ring R(m). By default the parent object S will depend only on R and m and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.The ResidueField constructor does the same thing as the ResidueRing constructor, but the resulting object has type belonging to Field rather than Ring, so it can be used anywhere a field is expected in AbstractAlgebra.jl. No check is made for maximality of the ideal generated by m.Here are some examples of creating residue rings and making use of the resulting parent objects to coerce various elements into the residue ring.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S()\ng = S(123)\nh = S(BigInt(1234))\nk = S(x + 1)All of the examples here are generic residue rings, but specialised implementations of residue rings provided by external modules will also usually provide a ResidueRing constructor to allow creation of their residue rings."
},

{
    "location": "residue/#Basic-ring-functionality-1",
    "page": "Generic residue rings",
    "title": "Basic ring functionality",
    "category": "section",
    "text": "Residue rings in AbstractAlgebra.jl implement the full Ring interface. Of course the entire Residue Ring interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nf = S(x + 1)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nm = modulus(S)\nU = base_ring(S)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)"
},

{
    "location": "residue/#Residue-ring-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic residue rings",
    "title": "Residue ring functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any residue ring module that implements the full Residue Ring interface. This includes AbstractAlgebra.jl\'s own generic residue rings.But if a C library provides all the functionality documented in the Residue Ring interface, then all the functions described here will also be automatically supplied by AbstractAlgebra.jl for that residue ring type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "residue/#AbstractAlgebra.Generic.modulus-Tuple{ResElem}",
    "page": "Generic residue rings",
    "title": "AbstractAlgebra.Generic.modulus",
    "category": "method",
    "text": "modulus(R::AbstractAlgebra.ResElem)\n\nReturn the modulus a of the residue ring S = R(a) that the supplied residue r belongs to.\n\n\n\n"
},

{
    "location": "residue/#AbstractAlgebra.Generic.isunit-Tuple{ResElem}",
    "page": "Generic residue rings",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "method",
    "text": "isunit(a::AbstractAlgebra.ResElem)\n\nReturn true if the supplied element a is invertible in the residue ring it belongs to, otherwise return false.\n\n\n\n"
},

{
    "location": "residue/#Basic-functionality-1",
    "page": "Generic residue rings",
    "title": "Basic functionality",
    "category": "section",
    "text": "modulus(::AbstractAlgebra.ResElem)isunit(::AbstractAlgebra.ResElem)ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R, x^3 + 3x + 1)\n\nr = S(x + 1)\n\na = modulus(S)\nisunit(r) == true"
},

{
    "location": "residue/#Base.inv-Tuple{ResElem}",
    "page": "Generic residue rings",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(a::AbstractAlgebra.ResElem)\n\nReturn the inverse of the element a in the residue ring. If an impossible inverse is encountered, an exception is raised.\n\n\n\n"
},

{
    "location": "residue/#Inversion-1",
    "page": "Generic residue rings",
    "title": "Inversion",
    "category": "section",
    "text": "inv(::AbstractAlgebra.ResElem)ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R)\n\nf = S(x + 1)\n\ng = inv(f)"
},

{
    "location": "residue/#Base.gcd-Union{Tuple{T}, Tuple{ResElem{T},ResElem{T}}} where T<:RingElem",
    "page": "Generic residue rings",
    "title": "Base.gcd",
    "category": "method",
    "text": "gcd(a::AbstractAlgebra.ResElem{T}, b::AbstractAlgebra.ResElem{T}) where {T <: RingElement}\n\nReturn a greatest common divisor of a and b if one exists. This is done by taking the greatest common divisor of the data associated with the supplied residues and taking its greatest common divisor with the modulus.\n\n\n\n"
},

{
    "location": "residue/#Greatest-common-divisor-1",
    "page": "Generic residue rings",
    "title": "Greatest common divisor",
    "category": "section",
    "text": "gcd{T <: RingElem}(::ResElem{T}, ::ResElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = ResidueRing(R)\n\nf = S(x + 1)\ng = S(x^2 + 2x + 1)\n\nh = gcd(f, g)"
},

{
    "location": "fields/#",
    "page": "Field Interface",
    "title": "Field Interface",
    "category": "page",
    "text": ""
},

{
    "location": "fields/#Field-Interface-1",
    "page": "Field Interface",
    "title": "Field Interface",
    "category": "section",
    "text": "AbstractAlgebra.jl generic code makes use of a standardised set of functions which it expects to be implemented for all fields. Here we document this interface. All libraries which want to make use of the generic capabilities of AbstractAlgebra.jl must supply all of the required functionality for their fields."
},

{
    "location": "fields/#Types-1",
    "page": "Field Interface",
    "title": "Types",
    "category": "section",
    "text": "Most fields must supply two types:a type for the parent object (representing the field itself)\na type for elements of that fieldFor example, the generic fraction field type in AbstractAlgebra.jl provides two  types in generic/GenericTypes.jl: Generic.FracField{T} for the parent objects\nGeneric.Frac{T} for the actual fractionsThe parent type must belong to AbstractAlgebra.Field and the element type must belong to AbstractAlgebra.FieldElem. Of course, the types may belong to these abstract types transitively.For parameterised fields, we advise that the types of both the parent objects and element objects to be parameterised by the types of the elements of the base ring.There can be variations on this theme: e.g. in some areas of mathematics there is a notion of a coefficient domain, in which case it may make sense to parameterise all types by the type of elements of this coefficient domain. But note that this may have implications for the ad hoc operators one might like to explicitly implement."
},

{
    "location": "fields/#Parent-object-caches-1",
    "page": "Field Interface",
    "title": "Parent object caches",
    "category": "section",
    "text": "In many cases, it is desirable to have only one object in the system to represent each field. This means that if the same field is constructed twice, elements of the two fields will be compatible as far as arithmetic is concerned.In order to facilitate this, global caches of fields are stored in AbstractAlgebra.jl, usually implemented using dictionaries. For example, the Generic.FracField parent objects are looked up in a dictionary FracDict to see if they have been previously defined.Whether these global caches are provided or not, depends on both mathematical and algorithmic considerations. E.g. in the case of number fields, it isn\'t desirable to identify all number fields with the same defining polynomial, as they may be considered with distinct embeddings into one another. In other cases, identifying whether two fields are the same may be prohibitively expensive. Generally, it may only make sense algorithmically to identify two fields if they were constructed from identical data.If a global cache is provided, it must be optionally possible to construct the parent objects without caching. This is done by passing a boolean value cached to the inner constructor of the parent object. See generic/GenericTypes.jl` for examples of how to construct and handle such caches."
},

{
    "location": "fields/#Required-functions-for-all-fields-1",
    "page": "Field Interface",
    "title": "Required functions for all fields",
    "category": "section",
    "text": "In the following, we list all the functions that are required to be provided for fields in AbstractAlgebra.jl or by external libraries wanting to use AbstractAlgebra.jl.We give this interface for fictitious types MyParent for the type of the field parent object R and MyElem for the type of the elements of the field.Note that generic functions in AbstractAlgebra.jl may not rely on the existence of functions that are not documented here. If they do, those functions will only be available for fields that implement that additional functionality, and should be documented as such.In the first place, all fields are rings and therefore any field type must implement all of the Ring interface. The functionality below is in addition to this basic functionality."
},

{
    "location": "fields/#Data-type-and-parent-object-methods-1",
    "page": "Field Interface",
    "title": "Data type and parent object methods",
    "category": "section",
    "text": "characteristic(R::MyParent)Return the characteristic of the field."
},

{
    "location": "fields/#Basic-manipulation-of-rings-and-elements-1",
    "page": "Field Interface",
    "title": "Basic manipulation of rings and elements",
    "category": "section",
    "text": "isunit(f::MyElem)Return true if the given element is invertible, i.e. nonzero in the field."
},

{
    "location": "fields/#Inversion-1",
    "page": "Field Interface",
    "title": "Inversion",
    "category": "section",
    "text": "inv(f::MyElem)Return the inverse of the given element in the field. If f = 0, an error is thrown."
},

{
    "location": "fraction_fields/#",
    "page": "Fraction Field Interface",
    "title": "Fraction Field Interface",
    "category": "page",
    "text": ""
},

{
    "location": "fraction_fields/#Fraction-Field-Interface-1",
    "page": "Fraction Field Interface",
    "title": "Fraction Field Interface",
    "category": "section",
    "text": "Fraction fields are supported in AbstractAlgebra.jl, at least for gcd domains. In addition to the standard Ring interface, some additional functions are required to be present for fraction fields."
},

{
    "location": "fraction_fields/#Types-and-parents-1",
    "page": "Fraction Field Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for fraction fields and their elements:FracField{T} is the abstract type for fraction field parent types\nFracElem{T} is the abstract type for types of fractionsWe have that FracField{T} <: AbstractAlgebra.Field and  FracElem{T} <: AbstractAlgebra.FieldElem.Note that both abstract types are parameterised. The type T should usually be the type of elements of the base ring of the fraction field.Fraction fields should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Fraction fields should at least be distinguished based on their base ring.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "fraction_fields/#Required-functionality-for-fraction-fields-1",
    "page": "Fraction Field Interface",
    "title": "Required functionality for fraction fields",
    "category": "section",
    "text": "In addition to the required functionality for the Field interface the Fraction Field interface has the following required functions.We suppose that R is a fictitious base ring, and that S is the fraction field with  parent object S of type MyFracField{T}. We also assume the fractions in the field  have type MyFrac{T}, where T is the type of elements of the base ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem."
},

{
    "location": "fraction_fields/#Constructors-1",
    "page": "Fraction Field Interface",
    "title": "Constructors",
    "category": "section",
    "text": "We provide the following constructors. Note that these constructors don\'t require construction of the parent object first. This is easier to achieve if the fraction element type doesn\'t contain a reference to the parent object, but merely contains a reference to the base ring. The parent object can then be constructed on demand.//(x::T, y::T) where T <: AbstractAlgebra.RingElemReturn the fraction xy.//(x::T, y::AbstractAlgebra.FracElem{T}) where T <: AbstractAlgebra.RingElemReturn xy where x is in the base ring of y.//(x::AbstractAlgebra.FracElem{T}, y::T) where T <: AbstractAlgebra.RingElemReturn xy where y is in the base ring of x.ExamplesR, x = PolynomialRing(ZZ, \"x\")\n\nf = (x^2 + x + 1)//(x^3 + 3x + 1)\ng = f//x\nh = x//f"
},

{
    "location": "fraction_fields/#Basic-manipulation-of-fields-and-elements-1",
    "page": "Fraction Field Interface",
    "title": "Basic manipulation of fields and elements",
    "category": "section",
    "text": "numerator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElemGiven a fraction d = ab return a, where ab is in lowest terms with respect to the canonical_unit and gcd functions on the base ring.denominator(d::MyFrac{T}) where T <: AbstractAlgebra.RingElemGiven a fraction d = ab return b, where ab is in lowest terms with respect to the canonical_unit and gcd functions on the base ring.ExamplesR, x = PolynomialRing(QQ, \"x\")\n\nf = (x^2 + x + 1)//(x^3 + 3x + 1)\n\nn = numerator(f)\nd = denominator(f)"
},

{
    "location": "fraction/#",
    "page": "Generic fraction fields",
    "title": "Generic fraction fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "fraction/#Generic-fraction-fields-1",
    "page": "Generic fraction fields",
    "title": "Generic fraction fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/generic/Fraction.jl for generic fraction fields over any gcd domain belonging to the AbstractAlgebra.jl abstract type hierarchy.As well as implementing the Fraction Field interface a number of generic algorithms are implemented for fraction fields. We describe this generic functionality below.All of the generic functionality is part of a submodule of AbstractAlgebra called Generic. This is exported by default so that it is not necessary to qualify the function names with the submodule name."
},

{
    "location": "fraction/#Types-and-parent-objects-1",
    "page": "Generic fraction fields",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Fractions implemented using the AbstractAlgebra generics have type Generic.Frac{T} where T is the type of elements of the base ring. See the file src/generic/GenericTypes.jl for details.Parent objects of such fraction elements have type Generic.FracField{T}.The fraction element types belong to the abstract type AbstractAlgebra.FracElem{T} and the fraction field types belong to the abstract type AbstractAlgebra.FracRing{T}. This enables one to write generic functions that can accept any AbstractAlgebra fraction type.Note that both the generic fraction field type Generic.FracField{T} and the abstract type it belongs to, AbstractAlgebra.FracField{T} are both called FracField. The  former is a (parameterised) concrete type for a fraction field over a given base ring whose elements have type T. The latter is an abstract type representing all fraction field types in AbstractAlgebra.jl, whether generic or very specialised (e.g. supplied by a C library)."
},

{
    "location": "fraction/#Fraction-field-constructors-1",
    "page": "Generic fraction fields",
    "title": "Fraction field constructors",
    "category": "section",
    "text": "In order to construct fractions in AbstractAlgebra.jl, one can first construct the fraction field itself. This is accomplished with the following constructor.FractionField(R::AbstractAlgebra.Ring; cached::Bool = true)Given a base ring R return the parent object of the fraction field of R. By default the parent object S will depend only on R and will be cached. Setting the optional argument cached to false will prevent the parent object S from being cached.Here are some examples of creating fraction fields and making use of the resulting parent objects to coerce various elements into the fraction field.ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS = FractionField(R)\n\nf = S()\ng = S(123)\nh = S(BigInt(1234))\nk = S(x + 1)All of the examples here are generic fraction fields, but specialised implementations of fraction fields provided by external modules will also usually provide a FractionField constructor to allow creation of the fraction fields they provide."
},

{
    "location": "fraction/#Basic-field-functionality-1",
    "page": "Generic fraction fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "Fraction fields in AbstractAlgebra.jl implement the full Field interface. Of course the entire Fraction Field interface is also implemented.We give some examples of such functionality.ExamplesR, x = PolynomialRing(QQ, \"x\")\nS = FractionField(R)\n\nf = S(x + 1)\ng = (x^2 + x + 1)//(x^3 + 3x + 1)\n\nh = zero(S)\nk = one(S)\nisone(k) == true\niszero(f) == false\nm = characteristic(S)\nU = base_ring(S)\nV = base_ring(f)\nT = parent(f)\nr = deepcopy(f)\nn = numerator(g)\nd = denominator(g)"
},

{
    "location": "fraction/#Fraction-field-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic fraction fields",
    "title": "Fraction field functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality listed below is automatically provided by AbstractAlgebra.jl for any fraction field module that implements the full Fraction Field interface. This includes AbstractAlgebra.jl\'s own generic fraction fields.But if a C library provides all the functionality documented in the Fraction Field interface, then all the functions described here will also be automatically supplied by AbstractAlgebra.jl for that fraction field type.Of course, modules are free to provide specific implementations of the functions described here, that override the generic implementation."
},

{
    "location": "fraction/#Base.gcd-Union{Tuple{T}, Tuple{FracElem{T},FracElem{T}}} where T<:RingElem",
    "page": "Generic fraction fields",
    "title": "Base.gcd",
    "category": "method",
    "text": "gcd(a::AbstractAlgebra.FracElem{T}, b::AbstractAlgebra.FracElem{T}) where {T <: RingElem}\n\nReturn a greatest common divisor of a and b if one exists. N.B: we define the GCD of ab and cd to be gcd(ad bc)bd, reduced to lowest terms. This requires the existence of a greatest common divisor function for the base ring.\n\n\n\n"
},

{
    "location": "fraction/#Greatest-common-divisor-1",
    "page": "Generic fraction fields",
    "title": "Greatest common divisor",
    "category": "section",
    "text": "gcd{T <: RingElem}(::FracElem{T}, ::FracElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\n\nf = (x + 1)//(x^3 + 3x + 1)\ng = (x^2 + 2x + 1)//(x^2 + x + 1)\n\nh = gcd(f, g)"
},

{
    "location": "fraction/#AbstractAlgebra.Generic.remove-Union{Tuple{T}, Tuple{FracElem{T},T}} where T<:RingElem",
    "page": "Generic fraction fields",
    "title": "AbstractAlgebra.Generic.remove",
    "category": "method",
    "text": "remove(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}\n\nReturn the tuple n x such that z = p^nx where x has valuation 0 at p.\n\n\n\n"
},

{
    "location": "fraction/#AbstractAlgebra.Generic.valuation-Union{Tuple{T}, Tuple{FracElem{T},T}} where T<:RingElem",
    "page": "Generic fraction fields",
    "title": "AbstractAlgebra.Generic.valuation",
    "category": "method",
    "text": "valuation(z::AbstractAlgebra.FracElem{T}, p::T) where {T <: RingElem}\n\nReturn the valuation of z at p.\n\n\n\n"
},

{
    "location": "fraction/#Remove-and-valuation-1",
    "page": "Generic fraction fields",
    "title": "Remove and valuation",
    "category": "section",
    "text": "When working over a Euclidean domain, it is convenient to extend valuations to the fraction field. To facilitate this, we define the following functions.remove{T <: RingElem}(::FracElem{T}, ::T)valuation{T <: RingElem}(::FracElem{T}, ::T)ExamplesR, x = PolynomialRing(ZZ, \"x\")\n\nf = (x + 1)//(x^3 + 3x + 1)\ng = (x^2 + 1)//(x^2 + x + 1)\n\nv, q = remove(f^3*g, x + 1)\nv = valuation(f^3*g, x + 1)"
},

{
    "location": "rational/#",
    "page": "Rational field",
    "title": "Rational field",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "rational/#Rational-field-1",
    "page": "Rational field",
    "title": "Rational field",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/Rational.jl for making Julia Rational{BigInt}s conform to the AbstractAlgebra.jl Field interface.In addition to providing a parent object QQ for Julia Rational{BigInt}s, we implement any additional functionality required by AbstractAlgebra.jl.Because Rational{BigInt} cannot be directly included in the AbstractAlgebra.jl abstract type hierarchy, we achieve integration of Julia Rational{BigInt}s by introducing a type union, called FieldElement, which is a union of AbstractAlgebra.FieldElem and a number of Julia types, including Rational{BigInt}. Everywhere that FieldElem is notionally used in AbstractAlgebra.jl, we are in fact using FieldElement, with additional care being taken to avoid ambiguities.The details of how this is done are technical, and we refer the reader to the implementation for details. For most intents and purposes, one can think of the Julia Rational{BigInt} type as belonging to AbstractAlgebra.FieldElem.One other technicality is that Julia defines certain functions for Rational{BigInt}, such as sqrt and exp differently to what AbstractAlgebra.jl requires. To get around this, we redefine these functions internally to AbstractAlgebra.jl, without redefining them for users of AbstractAlgebra.jl. This allows the internals of AbstractAlgebra.jl to function correctly, without broadcasting pirate definitions of already defined Julia functions to the world.To access the internal definitions, one can use AbstractAlgebra.sqrt and AbstractAlgebra.exp, etc."
},

{
    "location": "rational/#Types-and-parent-objects-1",
    "page": "Rational field",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Rationals have type Rational{BigInt}, as in Julia itself. We simply supplement the functionality for this type as required for computer algebra.The parent objects of such integers has type Rationals{BigInt}.For convenience, we also make Rational{Int} a part of the AbstractAlgebra.jl type hierarchy and its parent object (accessible as qq) has type Rationals{Int}. But we caution that this type is not particularly useful as a model of the rationals and may not function as expected within AbstractAlgebra.jl."
},

{
    "location": "rational/#Rational-constructors-1",
    "page": "Rational field",
    "title": "Rational constructors",
    "category": "section",
    "text": "In order to construct rationals in AbstractAlgebra.jl, one can first construct the rational field itself. This is accomplished using either of the following constructors.FractionField(R::Integers{BigInt})Rationals{BigInt}()This gives the unique object of type Rationals{BigInt} representing the field of rationals in AbstractAlgebra.jl.In practice, one simply uses QQ which is assigned to be the return value of the above constructor. There is no need to call the constructor in practice.Here are some examples of creating the rational field and making use of the resulting parent object to coerce various elements into the field.Examplesf = QQ()\ng = QQ(123)\nh = QQ(BigInt(1234))\nk = QQ(BigInt(12), BigInt(7))\n\nQQ == FractionField(ZZ)"
},

{
    "location": "rational/#Basic-field-functionality-1",
    "page": "Rational field",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The rational field in AbstractAlgebra.jl implements the full Field and Fraction Field interfaces.We give some examples of such functionality.Examplesf = QQ(12, 7)\n\nh = zero(QQ)\nk = one(QQ)\nisone(k) == true\niszero(f) == false\nU = base_ring(QQ)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)\ng = f + 12\nr = ZZ(12)//ZZ(7)\nn = numerator(r)"
},

{
    "location": "rational/#Rational-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Rational field",
    "title": "Rational functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "The functionality below supplements that provided by Julia itself for its Rational{BigInt} type."
},

{
    "location": "rational/#AbstractAlgebra.sqrt-Tuple{Rational{BigInt}}",
    "page": "Rational field",
    "title": "AbstractAlgebra.sqrt",
    "category": "method",
    "text": "sqrt(a::Rational{T}) where T <: Integer\n\nReturn the square root of a if it is the square of a rational, otherwise throw an error.\n\n\n\n"
},

{
    "location": "rational/#AbstractAlgebra.exp-Tuple{Rational{BigInt}}",
    "page": "Rational field",
    "title": "AbstractAlgebra.exp",
    "category": "method",
    "text": "exp(a::Rational{T}) where T <: Integer\n\nReturn 1 if a = 0, otherwise throw an exception.\n\n\n\n"
},

{
    "location": "rational/#Square-root-1",
    "page": "Rational field",
    "title": "Square root",
    "category": "section",
    "text": "AbstractAlgebra.sqrt(a::Rational{BigInt})AbstractAlgebra.exp(a::Rational{BigInt})Examplesd = AbstractAlgebra.sqrt(ZZ(36)//ZZ(25))\nm = AbstractAlgebra.exp(ZZ(0)//ZZ(1))"
},

{
    "location": "finfield/#",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "finfield/#Finite-fields-1",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/GF.jl for finite fields. The module is a naive implementation that supports only fields of degree 1 (prime fields). They are modelled as mathbbZpmathbbZ for p a prime."
},

{
    "location": "finfield/#Types-and-parent-objects-1",
    "page": "Finite fields",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Finite fields have type GFField{T} where T is either Int or BigInt.Elements of such a finite field have type gfelem{T}."
},

{
    "location": "finfield/#AbstractAlgebra.GF-Union{Tuple{T}, Tuple{T}} where T<:Integer",
    "page": "Finite fields",
    "title": "AbstractAlgebra.GF",
    "category": "method",
    "text": "GF(p::T) where T <: Integer\n\nReturn the finite field mathbbF_p, where p is a prime. The integer p is not checked for primality, but the behaviour of the resulting object is undefined if p is composite.\n\n\n\n"
},

{
    "location": "finfield/#Finite-field-constructors-1",
    "page": "Finite fields",
    "title": "Finite field constructors",
    "category": "section",
    "text": "In order to construct finite fields in AbstractAlgebra.jl, one must first construct the field itself. This is accomplished with the following constructors.GF(p::T) where T <: IntegerHere are some examples of creating a finite field and making use of the resulting parent object to coerce various elements into the field.ExamplesF = GF(13)\n\ng = F(3)\nh = F(g)"
},

{
    "location": "finfield/#Basic-field-functionality-1",
    "page": "Finite fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The finite field module in AbstractAlgebra.jl implements the full Field interface.We give some examples of such functionality.ExamplesF = GF(13)\n\nh = zero(F)\nk = one(F)\nisone(k) == true\niszero(f) == false\nU = base_ring(F)\nV = base_ring(h)\nT = parent(h)\nh == deepcopy(h)\nh = h + 2\nm = inv(k)"
},

{
    "location": "finfield/#AbstractAlgebra.Generic.gen-Union{Tuple{GFField{T}}, Tuple{T}} where T<:Integer",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.gen",
    "category": "method",
    "text": "gen(R::GFField{T}) where T <: Integer\n\nReturn a generator of the field. Currently this returns 1.\n\n\n\n"
},

{
    "location": "finfield/#AbstractAlgebra.Generic.order-Tuple{AbstractAlgebra.GFField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.order",
    "category": "method",
    "text": "order(R::GFField)\n\nReturn the order, i.e. the number of element in, the given finite field.\n\n\n\n"
},

{
    "location": "finfield/#AbstractAlgebra.Generic.degree-Tuple{AbstractAlgebra.GFField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.degree",
    "category": "method",
    "text": "degree(R::GFField)\n\nReturn the degree of the given finite field.\n\n\n\n"
},

{
    "location": "finfield/#Basic-manipulation-of-fields-and-elements-1",
    "page": "Finite fields",
    "title": "Basic manipulation of fields and elements",
    "category": "section",
    "text": "gen{T <: Integer}(F::GFField{T})order(F::GFField)degree(F::GFField)ExamplesF = GF(13)\n\nd = degree(F)\nn = order(F)\ng = gen(F)"
},

{
    "location": "real/#",
    "page": "Real field",
    "title": "Real field",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "real/#Real-field-1",
    "page": "Real field",
    "title": "Real field",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a module, implemented in src/julia/Float.jl for making Julia BigFloats conform to the AbstractAlgebra.jl Field interface.In addition to providing a parent object RealField for Julia BigFloats, we implement any additional functionality required by AbstractAlgebra.jl.Because BigFloat cannot be directly included in the AbstractAlgebra.jl abstract type hierarchy, we achieve integration of Julia BigFloats by introducing a type union, called FieldElement, which is a union of AbstractAlgebra.FieldElem and a number of Julia types, including BigFloat. Everywhere that FieldElem is notionally used in AbstractAlgebra.jl, we are in fact using FieldElement, with additional care being taken to avoid ambiguities.The details of how this is done are technical, and we refer the reader to the implementation for details. For most intents and purposes, one can think of the Julia BigFloat type as belonging to AbstractAlgebra.FieldElem."
},

{
    "location": "real/#Types-and-parent-objects-1",
    "page": "Real field",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Reals have type BigFloat, as in Julia itself. We simply supplement the functionality for this type as required for computer algebra.The parent objects of such integers has type Floats{BigFloat}.For convenience, we also make Float64 a part of the AbstractAlgebra.jl type hierarchy and its parent object (accessible as RDF) has type Floats{Float64}."
},

{
    "location": "real/#Rational-constructors-1",
    "page": "Real field",
    "title": "Rational constructors",
    "category": "section",
    "text": "In order to construct reals in AbstractAlgebra.jl, one can first construct the real field itself. This is accomplished using the following constructor.Floats{BigFloat}()This gives the unique object of type Floats{BigFloat} representing the field of reals in AbstractAlgebra.jl.In practice, one simply uses RealField which is assigned to be the return value of the above constructor. There is no need to call the constructor in practice.Here are some examples of creating the real field and making use of the resulting parent object to coerce various elements into the field.ExamplesRR = RealField\n\nf = RR()\ng = RR(123)\nh = RR(BigInt(1234))\nk = RR(12//7)\nm = RR(2.3)"
},

{
    "location": "real/#Basic-field-functionality-1",
    "page": "Real field",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The real field in AbstractAlgebra.jl implements the full Field interface.We give some examples of such functionality.ExamplesRR = RealField\nf = RR(12//7)\n\nh = zero(RR)\nk = one(RR)\nisone(k) == true\niszero(f) == false\nU = base_ring(RR)\nV = base_ring(f)\nT = parent(f)\nf == deepcopy(f)\ng = f + 12\nm = inv(g)"
},

{
    "location": "numberfield/#",
    "page": "Number fields",
    "title": "Number fields",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "numberfield/#Number-fields-1",
    "page": "Number fields",
    "title": "Number fields",
    "category": "section",
    "text": "AbstractAlgebra.jl provides a very naive implementation of number fields. This allows arithmetic in algebraic number fields, which are currently modeled as mathbbQx modulo an irreducible polynomial, i.e. as a residue field.In fact, the definition of the number field constructor is currently given in src/generic/ResidueField.jl and no type is defined for a number field. The definition mainly exists for testing purposes. It may later be replaced by a more standard implementation. For a more fully fleshed out number field implementation (based on a very high performance C library), see Nemo.jl."
},

{
    "location": "numberfield/#Number-field-constructors-1",
    "page": "Number fields",
    "title": "Number field constructors",
    "category": "section",
    "text": "In order to construct number fields in AbstractAlgebra.jl, one must first construct the field itself. This is accomplished with the following constructor.NumberField(f::AbstractAlgebra.Generic.Poly{Rational{BigInt}}, s::AbstractString, t = \"\\$\"; cached = true)Given an irreducible defining polynomial f in mathbbQx, return a tuple (K x) consisting of the number field defined by that polynomial and a generator. The string fields are currently ignored, but are reserved for future use.Currently the generator of the number field prints the same way as the variable in mathbbQx.ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\n\nf = a^2 + 2a + 7"
},

{
    "location": "numberfield/#Basic-field-functionality-1",
    "page": "Number fields",
    "title": "Basic field functionality",
    "category": "section",
    "text": "The number field module in AbstractAlgebra.jl implements the full Field and ResidueRing interfaces."
},

{
    "location": "perm/#",
    "page": "Permutations and Permutation groups",
    "title": "Permutations and Permutation groups",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra\nDocTestSetup = quote\n    using AbstractAlgebra\nend"
},

{
    "location": "perm/#AbstractAlgebra.Generic.setpermstyle",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.setpermstyle",
    "category": "function",
    "text": "setpermstyle(format::Symbol)\n\nSelect the style in which permutations are displayed (in REPL or in general as string). This can be either:array - as vectors of integers whose n-th position represents thevalue at n), or:cycles - as, more familiar for mathematicians, decomposition intodisjoint cycles, where the value at n is represented by the entry immediately following n in a cycle (the default).\n\nThe difference is purely esthetical.\n\nExamples:\n\njulia> Generic.setpermstyle(:array)\n:array\n\njulia> perm([2,3,1,5,4])\n[2, 3, 1, 5, 4]\n\njulia> Generic.setpermstyle(:cycles)\n:cycles\n\njulia> perm([2,3,1,5,4])\n(1,2,3)(4,5)\n\n\n\n"
},

{
    "location": "perm/#Permutations-and-Permutation-groups-1",
    "page": "Permutations and Permutation groups",
    "title": "Permutations and Permutation groups",
    "category": "section",
    "text": "AbstractAlgebra.jl provides rudimentary native support for permutation groups (implemented in src/generic/PermGroups.jl). All functionality of permutations is accesible in the Generic submodule.Permutations are represented internally via vector of integers, wrapped in type perm{T}, where T<:Integer carries the information on the type of elements of a permutation. Permutation groups are singleton parent objects of type PermGroup{T} and are used mostly to store the length of a permutation, since it is not included in the permutation type.Permutation groups are created using the PermGroup (inner) constructor. However, for convenience we definePermutationGroup = PermGroupso that permutation groups can be created using PermutationGroup instead of PermGroup.Both PermGroup and perm and can be parametrized by any type T<:Integer . By default the parameter is the Int-type native to the systems architecture. However, if you are sure that your permutations are small enough to fit into smaller integer type (such as Int32, Uint16, or even Int8), you may choose to change the parametrizing type accordingly. In practice this may result in decreased memory footprint (when storing multiple permutations) and noticable faster performance, if your workload is heavy in operations on permutations, which e.g. does not fit into cache of your cpu.All the permutation group types belong to the Group abstract type and the corresponding permutation element types belong to the GroupElem abstract type.Generic.setpermstyle"
},

{
    "location": "perm/#AbstractAlgebra.Generic.perm",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.perm",
    "category": "type",
    "text": "perm{T<:Integer}\n\nThe type of permutations. Fieldnames:d::Vector{T} - vector representing the permutation\nmodified::Bool - bit to check the validity of cycle decomposition\ncycles::CycleDec{T} - (cached) cycle decompositionPermutation p consists of a vector (p.d) of n integers from 1 to n. If the i-th entry of the vector is j, this corresponds to p sending i to j. The cycle decomposition (p.cycles) is computed on demand and should never be accessed directly. Use cycles(p) instead.There are two inner constructors of perm:perm(n::T) constructs the trivial perm{T}-permutation of length n.\nperm(v::Vector{T<:Integer}[,check=true]) constructs a permutationrepresented by v. By default perm constructor checks if the vector constitutes a valid permutation. To skip the check call perm(v, false).\n\nExamples:\n\njulia> perm([1,2,3])\n()\n\njulia> g = perm(Int32[2,3,1])\n(1,2,3)\n\njulia> typeof(g)\nAbstractAlgebra.Generic.perm{Int32}\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.PermGroup",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.PermGroup",
    "category": "type",
    "text": "PermGroup{T<:Integer}\n\nThe permutation group singleton type. PermGroup(n) constructs the permutation group S_n on n-symbols. The type of elements of the group is inferred from the type of n.\n\nExamples:\n\njulia> G = PermGroup(5)\nPermutation group over 5 elements\n\njulia> elem_type(G)\nAbstractAlgebra.Generic.perm{Int64}\n\njulia> H = PermGroup(UInt16(5))\nPermutation group over 5 elements\n\njulia> elem_type(H)\nAbstractAlgebra.Generic.perm{UInt16}\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.@perm_str",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.@perm_str",
    "category": "macro",
    "text": "perm\"...\"\n\nString macro to parse disjoint cycles into perm{Int}.Strings for the output of GAP could be copied directly into perm\"...\". Cycles of length 1 are not necessary, but could be included. A permutation of the minimal support is constructed, i.e. the maximal n in the decomposition determines the parent group S_n.\n\nExamples:\n\njulia> p = perm\"(1,3)(2,4)\"\n(1,3)(2,4)\n\njulia> typeof(p)\nAbstractAlgebra.Generic.perm{Int64}\n\njulia> parent(p) == PermutationGroup(4)\ntrue\n\njulia> p = perm\"(1,3)(2,4)(10)\"\n(1,3)(2,4)\n\njulia> parent(p) == PermutationGroup(10)\ntrue\n\n\n\n"
},

{
    "location": "perm/#Permutations-constructors-1",
    "page": "Permutations and Permutation groups",
    "title": "Permutations constructors",
    "category": "section",
    "text": "There are several methods to to construct permutations in AbstractAlgebra.jl.The easiest way is to directly call to the perm (inner) constructor:Generic.permSince the parent object can be reconstructed from the permutation itself, you can work with permutations without explicitely constructing the parent object.The other way is to first construct the permutation group they belong to. This is accomplished with the inner constructor PermGroup(n::Integer) which constructs the permutation group on n symbols and returns the parent object representing the group.Generic.PermGroupA vector of integers can be then coerced to a permutation via call to parent. The advantage is that the vector is automatically converted to the integer type fixed at the creation of the parent object.Examples:julia> G = PermutationGroup(BigInt(5)); p = G([2,3,1,5,4])\n(1,2,3)(4,5)\n\njulia> typeof(p)\nAbstractAlgebra.Generic.perm{BigInt}\n\njulia> H = PermutationGroup(UInt16(5)); r = H([2,3,1,5,4])\n(1,2,3)(4,5)\n\njulia> typeof(r)\nAbstractAlgebra.Generic.perm{UInt16}\n\njulia> H()\n()By default the coercion checks for non-unique values in the vector, but this can be switched off with G([2,3,1,5,4], false).Finally there is a perm\"...\" string macro to construct permutation from string input.@perm_str"
},

{
    "location": "perm/#Base.parent-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "Base.parent",
    "category": "method",
    "text": "parent(g::perm{T}) where T = PermGroup\n\nReturn the parent of the permutation g.\n\njulia> G = PermutationGroup(5); g = perm([3,4,5,2,1])\n(1,3,5)(2,4)\n\njulia> parent(g) == G\ntrue\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.elem_type-Tuple{AbstractAlgebra.PermGroup}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.elem_type",
    "category": "method",
    "text": "elem_type(::Type{PermGroup{T}}) where T = perm{T}\n\nReturn the type of elements of a permutation group.\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.parent_type-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.parent_type",
    "category": "method",
    "text": "parent_type(::Type{perm{T}}) where T = PermGroup{T}\n\nReturn the type of the parent of a permutation.\n\n\n\n"
},

{
    "location": "perm/#Permutation-interface-1",
    "page": "Permutations and Permutation groups",
    "title": "Permutation interface",
    "category": "section",
    "text": "The following basic functionality is provided by the default permutation group implementation in AbstractAlgebra.jl, to support construction of other generic constructions over permutation groups. Any custom permutation group implementation in AbstractAlgebra.jl should provide these functions along with the usual group element arithmetic and comparison.parent(::perm)\nelem_type(::PermGroup)\nparent_type(::perm)A custom implementation also needs to implement hash(::perm, ::UInt) and (possibly) deepcopy_internal(::perm, ::ObjectIdDict).note: Note\nPermutation group elements are mutable and so returning shallow copies is not sufficient.getindex(a::perm, n::Int)Allows access to entry n of the given permutation via the syntax a[n]. Note that entries are 1-indexed.setindex!(a::perm, d::Int, n::Int)Set the n-th entry of the given permutation to d. This allows Julia to provide the syntax a[n] = d for setting entries of a permutation. Entries are 1-indexed.note: Note\nUsing setindex! invalidates cycle decomposition cached in a permutation, i.e. it will be computed the next time cycle decomposition is needed.Given the parent object G for a permutation group, the following coercion functions are provided to coerce various arguments into the permutation group. Developers provide these by overloading the permutation group parent objects.G()Return the identity permutation.G(A::Vector{<:Integer})Return the permutation whose entries are given by the elements of the supplied vector.G(p::perm)Take a permutation that is already in the permutation group and simply return it. A copy of the original is not made if not necessary."
},

{
    "location": "perm/#AbstractAlgebra.Generic.cycles-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.cycles",
    "category": "method",
    "text": "cycles(g::perm{T}) where T<:Integer\n\nDecompose permutation g into disjoint cycles.Returns a CycleDec object which iterates over disjoint cycles of g. The ordering of cycles is not guaranteed, and the order within each cycle is computed up to a cyclic permutation. The cycle decomposition is cached in g and used in future computation of permtype, parity, sign, order and ^ (powering).\n\nExamples:\n\njulia> g = perm([3,4,5,2,1,6])\n(1,3,5)(2,4)\n\njulia> collect(cycles(g))\n3-element Array{Array{Int64,1},1}:\n [1, 3, 5]\n [2, 4]\n [6]\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.parity-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.parity",
    "category": "method",
    "text": "parity(g::perm{T}) where T\n\nReturn the parity of the given permutation, i.e. the parity of the number of transpositions in any decomposition of g into transpositions.parity returns 1 if the number is odd and 0 otherwise. parity uses cycle decomposition of g if already available, but will not compute it on demand. Since cycle structure is cached in g you may call cycles(g) before calling parity.\n\nExamples:\n\njulia> g = perm([3,4,1,2,5])\n(1,3)(2,4)\n\njulia> parity(g)\n0\n\njulia> g = perm([3,4,5,2,1,6])\n(1,3,5)(2,4)\n\njulia> parity(g)\n1\n\n\n\n"
},

{
    "location": "perm/#Base.sign-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "Base.sign",
    "category": "method",
    "text": "sign(g::perm{T}) where T\n\nReturn the sign of permutation.sign returns 1 if g is even and -1 if g is odd. sign represents the homomorphism from the permutation group to the unit group of mathbbZ whose kernel is the alternating group.\n\nExamples:\n\njulia> g = perm([3,4,1,2,5])\n(1,3)(2,4)\n\njulia> sign(g)\n1\n\njulia> g = perm([3,4,5,2,1,6])\n(1,3,5)(2,4)\n\njulia> sign(g)\n-1\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.permtype-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.permtype",
    "category": "method",
    "text": "permtype(g::perm)\n\nReturn the type of permutation g, i.e. lengths of disjoint cycles in cycle decomposition of g.The lengths are sorted in decreasing order by default. permtype(g) fully determines the conjugacy class of g.\n\nExamples:\n\njulia> g = perm([3,4,5,2,1,6])\n(1,3,5)(2,4)\n\njulia> permtype(g)\n3-element Array{Int64,1}:\n 3\n 2\n 1\n\njulia> G = PermGroup(5); e = parent(g)()\n()\n\njulia> permtype(e)\n6-element Array{Int64,1}:\n 1\n 1\n 1\n 1\n 1\n 1\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.order-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.order",
    "category": "method",
    "text": "order(a::perm) -> BigInt\n\nReturn the order of permutation a as BigInt.If you are sure that computation over T (or its Int promotion) will not overflow you may use the method order(T::Type, a::perm) which bypasses computation with BigInts and returns promote(T, Int).\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.order-Tuple{AbstractAlgebra.Generic.PermGroup}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.order",
    "category": "method",
    "text": "order(G::PermGroup) -> BigInt\n\nReturn the order of the full permutation group as BigInt.\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.elements!-Tuple{AbstractAlgebra.Generic.PermGroup}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.elements!",
    "category": "method",
    "text": "Generic.elements!(G::PermGroup)\n\nReturn an unsafe iterator over all permutations in G. Only one permutation is allocated and then modified in-place using the non-recursive Heaps algorithm.Note: you need to explicitely copy permutations intended to be stored or modified.\n\nExamples:\n\njulia> elts = Generic.elements!(PermGroup(5));\n\njulia> length(elts)\n120\n\njulia> for p in Generic.elements!(PermGroup(3))\n         println(p)\n       end\n()\n(1,2)\n(1,3,2)\n(2,3)\n(1,2,3)\n(1,3)\n\njulia> A = collect(Generic.elements!(PermGroup(3))); A\n6-element Array{AbstractAlgebra.Generic.perm{Int64},1}:\n (1,3)\n (1,3)\n (1,3)\n (1,3)\n (1,3)\n (1,3)\n\njulia> unique(A)\n1-element Array{AbstractAlgebra.Generic.perm{Int64},1}:\n (1,3)\n\n\n\n"
},

{
    "location": "perm/#Basic-manipulation-1",
    "page": "Permutations and Permutation groups",
    "title": "Basic manipulation",
    "category": "section",
    "text": "Numerous functions are provided to manipulate permutation group elements.cycles(::perm)Cycle structure is cached in a permutation, since once available, it provides a convenient shortcut in many other algorithms.parity(::perm)\nsign(::perm)\npermtype(::perm)\norder(::perm)\norder(::Generic.PermGroup)Note that even an Int64 can be easily overflowed when computing with permutation groups. Thus, by default, order returns (always correct) BigInts. If you are sure that the computation will not overflow, you may use order(::Type{T}, ...) to perform computations with machine integers. Julias standard promotion rules apply for the returned value.Since PermGroup implements the iterator protocole You may iterate over all permutations via simplefor p in PermutationGroup(n)\n   ...\nendIteration over all permutations in reasonable time, (i.e. in terms of minutes) is possible when n  13.You may also use the non-allocating Generic.elements! function for n  14 (or even 15 if you are patient enough), which is an order of mangitude faster.Generic.elements!(::Generic.PermGroup)However, since all permutations yielded by elements! are aliased (modified \"in-place\"), collect(Generic.elements!(PermGroup(n))) returns a vector of identical permutations.note: Note\nIf you intend to use or store elements yielded by elements! you need to deepcopy them explicitely."
},

{
    "location": "perm/#Base.:*-Union{Tuple{T}, Tuple{perm{T},perm{T}}} where T",
    "page": "Permutations and Permutation groups",
    "title": "Base.:*",
    "category": "method",
    "text": "*(g::perm{T}, h::perm{T}) where T\n\nReturn the composition h  g of two permutations.This corresponds to the action of permutation group on the set [1..n] on the right and follows the convention of GAP.If g and h are parametrized by different types, the result is promoted accordingly.\n\nExamples:\n\njulia> perm([2,3,1,4])*perm([1,3,4,2]) # (1,2,3)*(2,3,4)\n(1,3)(2,4)\n\n\n\n"
},

{
    "location": "perm/#Base.:^-Tuple{perm,Integer}",
    "page": "Permutations and Permutation groups",
    "title": "Base.:^",
    "category": "method",
    "text": "^(g::perm{T}, n::Integer) where T\n\nReturn the n-th power of a permutation g.By default g^n is computed by cycle decomposition of g if n > 3. Generic.power_by_squaring provides a different method for powering which may or may not be faster, depending on the particuar case. Due to caching of the cycle structure, repeated powering of g will be faster with the default method.\n\nExamples:\n\njulia> g = perm([2,3,4,5,1])\n(1,2,3,4,5)\n\njulia> g^3\n(1,4,2,5,3)\n\njulia> g^5\n()\n\n\n\n"
},

{
    "location": "perm/#Base.inv-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(g::perm)\n\nReturn the inverse of the given permutation, i.e. the permuation g^-1 such that g  g^-1 = g^-1  g is the identity permutation.\n\n\n\n"
},

{
    "location": "perm/#Arithmetic-operators-1",
    "page": "Permutations and Permutation groups",
    "title": "Arithmetic operators",
    "category": "section",
    "text": "*(::perm{T}, ::perm{T}) where T\n^(::perm, n::Integer)\ninv(::perm)Permutations parametrized by different types can be multiplied, and follow the standard julia integer promotion rules:g = rand(PermGroup(Int8(5)));\nh = rand(PermGroup(UInt32(5)));\ntypeof(g*h)\n\n# output\nAbstractAlgebra.Generic.perm{Int64}"
},

{
    "location": "perm/#Coercion-1",
    "page": "Permutations and Permutation groups",
    "title": "Coercion",
    "category": "section",
    "text": "The following coercions are available for G::PermGroup parent objects. Each of the methods perform basic sanity checks on the input which can be switched off by the second argument.Examples(G::PermGroup)()Return the identity element of G.(G::PermGrup)(::Vector{<:Integer}[, check=true])Turn a vector od integers into a permutation (performing conversion, if necessary).(G::PermGroup)(::perm{<:Integer}[, check=true])Coerce a permutation p into group G (performing the conversion, if necessary). If p is already an element of G no copy is performed.(G::PermGroup)(::String[, check=true])Parse the string input e.g. copied from the output of GAP. The method uses the same logic as perm\"...\" macro. The string is sanitized and checked for disjoint cycles. Both string(p::perm) (if setpermstyle(:cycles)) and string(cycles(p::perm)) are valid input for this method.(G::PermGroup{T})(::CycleDec{T}[, check=true]) where TTurn a cycle decomposition object into a permutation."
},

{
    "location": "perm/#Base.:==-Tuple{perm,perm}",
    "page": "Permutations and Permutation groups",
    "title": "Base.:==",
    "category": "method",
    "text": "==(g::perm, h::perm)\n\nReturn true if permutations are equal, otherwise return false.Permutations parametrized by different integer types are considered equal if they define the same permutation in the abstract permutation group.\n\nExamples:\n\njulia> g = perm(Int8[2,3,1])\n(1,2,3)\n\njulia> h = perm\"(3,1,2)\"\n(1,2,3)\n\njulia> g == h\ntrue\n\n\n\n"
},

{
    "location": "perm/#Base.:==-Tuple{AbstractAlgebra.Generic.PermGroup,AbstractAlgebra.Generic.PermGroup}",
    "page": "Permutations and Permutation groups",
    "title": "Base.:==",
    "category": "method",
    "text": "==(G::PermGroup, H::PermGroup)\n\nReturn true if permutation groups are equal, otherwise return false.Permutation groups on the same number of letters, but parametrized by different integer types are considered different.\n\nExamples:\n\njulia> G = PermGroup(UInt(5))\nPermutation group over 5 elements\n\njulia> H = PermGroup(5)\nPermutation group over 5 elements\n\njulia> G == H\nfalse\n\n\n\n"
},

{
    "location": "perm/#Comparison-1",
    "page": "Permutations and Permutation groups",
    "title": "Comparison",
    "category": "section",
    "text": "==(::perm, ::perm)\n==(::Generic.PermGroup, ::Generic.PermGroup)"
},

{
    "location": "perm/#Base.rand-Tuple{AbstractAlgebra.Generic.PermGroup}",
    "page": "Permutations and Permutation groups",
    "title": "Base.rand",
    "category": "method",
    "text": "rand(G::PermGroup{T}) where {T}\n\nReturn a random permutation from G.\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.matrix_repr-Tuple{perm}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.matrix_repr",
    "category": "method",
    "text": "matrix_repr(a::perm{T}) where T\n\nReturn the permutation matrix as sparse matrix representing a via natural embedding of the permutation group into general linear group over mathbbZ.\n\nExamples:\n\njulia> p = perm([2,3,1])\n(1,2,3)\n\njulia> matrix_repr(p)\n3×3 SparseMatrixCSC{Int64,Int64} with 3 stored entries:\n  [3, 1]  =  1\n  [1, 2]  =  1\n  [2, 3]  =  1\n\njulia> full(ans)\n3×3 Array{Int64,2}:\n 0  1  0\n 0  0  1\n 1  0  0\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.emb-Tuple{AbstractAlgebra.Generic.PermGroup,Array{Int64,1},Bool}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.emb",
    "category": "method",
    "text": "emb(G::PermGroup, V::Vector{Int}, check::Bool=true)\n\nReturn the natural embedding of a permutation group into G as the subgroup permuting points indexed by V.\n\nExamples:\n\njulia> p = perm([2,3,1])\n(1,2,3)\n\njulia> f = Generic.emb(PermGroup(5), [3,2,5]);\n\njulia> f(p)\n(2,5,3)\n\n\n\n"
},

{
    "location": "perm/#AbstractAlgebra.Generic.emb!-Tuple{perm,perm,Any}",
    "page": "Permutations and Permutation groups",
    "title": "AbstractAlgebra.Generic.emb!",
    "category": "method",
    "text": "emb!(result::perm, p::perm, V)\n\nEmbed permutation p into permutation result on the indices given by V.This corresponds to the natural embedding of S_k into S_n as the subgroup permuting points indexed by V.\n\nExamples:\n\njulia> p = perm([2,1,4,3])\n(1,2)(3,4)\n\njulia> Generic.emb!(perm(collect(1:5)), p, [3,1,4,5])\n(1,3)(4,5)\n\n\n\n"
},

{
    "location": "perm/#Misc-1",
    "page": "Permutations and Permutation groups",
    "title": "Misc",
    "category": "section",
    "text": "rand(::Generic.PermGroup)\nGeneric.matrix_repr(::perm)\nGeneric.emb(::Generic.PermGroup, ::Vector{Int}, ::Bool)\nGeneric.emb!(::perm, ::perm, V)"
},

{
    "location": "ytabs/#",
    "page": "Partitions and Young tableaux",
    "title": "Partitions and Young tableaux",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra\nDocTestSetup = quote\n    using AbstractAlgebra\nend\nDocTestFilters = r\"[0-9\\.]+ seconds \\(.*\\)\""
},

{
    "location": "ytabs/#Partitions-and-Young-tableaux-1",
    "page": "Partitions and Young tableaux",
    "title": "Partitions and Young tableaux",
    "category": "section",
    "text": "AbstractAlgebra.jl provides basic support for computations with Young tableaux, skew diagrams and the characters of permutation groups (implemented src/generic/YoungTabs.jl). All functionality of permutations is accesible in the Generic submodule."
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.Partition",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.Partition",
    "category": "type",
    "text": "Partition(part::Vector{<:Integer}[, check::Bool=true]) <: AbstractVector{Int}\n\nRepresent integer partition in the non-increasing order.part will be sorted, if necessary. Checks for validity of input can be skipped by calling the (inner) constructor with false as the second argument.Functionally Partition is a thin wrapper over Vector{Int}.\n\nFieldnames:n::Int - the partitioned number\npart::Vector{Int} - a non-increasing sequence of summands of n.\n\nExamples:\n\njulia> p = Partition([4,2,1,1,1])\n4₁2₁1₃\n\njulia> p.n == sum(p.part)\ntrue\n\n\n\n"
},

{
    "location": "ytabs/#Partitions-1",
    "page": "Partitions and Young tableaux",
    "title": "Partitions",
    "category": "section",
    "text": "The basic underlying object for those concepts is Partition of a number n, i.e. a sequence of positive integers n_1 ldots n_k which sum to n. Partitions in AbstractAlgebra.jl are represented internally by non-increasing Vectors of Ints. Partitions are printed using the standard notation, i.e. 9 = 4 + 2 + 1 + 1 + 1 is shown as 4_1 2_1 1_3 with the subscript indicating the count of a summand in the partition.Generic.Partition"
},

{
    "location": "ytabs/#Base.size-Tuple{AbstractAlgebra.Generic.Partition}",
    "page": "Partitions and Young tableaux",
    "title": "Base.size",
    "category": "method",
    "text": "size(p::Partition)\n\nReturn the size of the vector which represents the partition.\n\nExamples:\n\njulia> p = Partition([4,3,1]); size(p)\n(3,)\n\n\n\n"
},

{
    "location": "ytabs/#Base.getindex-Tuple{AbstractAlgebra.Generic.Partition,Integer}",
    "page": "Partitions and Young tableaux",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(p::Partition, i::Integer)\n\nReturn the i-th part (in decreasing order) of the partition.\n\n\n\n"
},

{
    "location": "ytabs/#Base.setindex!-Tuple{AbstractAlgebra.Generic.Partition,Integer,Integer}",
    "page": "Partitions and Young tableaux",
    "title": "Base.setindex!",
    "category": "method",
    "text": "setindex!(p::Partition, v::Integer, i::Integer)\n\nSet the i-th part of partition p to v. setindex! will throw an error if the operation violates the non-increasing assumption.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.AllParts",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.AllParts",
    "category": "type",
    "text": "AllParts(n::Integer)\n\nReturn an iterator over all integer Partitions of n. Partitions are produced in ascending order according to RuleAsc (Algorithm 3.1) fromJerome Kelleher and Barry O’Sullivan, Generating All Partitions: A Comparison Of Two Encodings ArXiv:0909.2331See also Combinatorics.partitions(1:n).\n\nExamples\n\njulia> ap = AllParts(5);\n\njulia> collect(ap)\n7-element Array{AbstractAlgebra.Generic.Partition,1}:\n 1₅\n 2₁1₃\n 3₁1₂\n 2₂1₁\n 4₁1₁\n 3₁2₁\n 5₁\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic._numpart",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic._numpart",
    "category": "function",
    "text": "_numpart(n::Integer)\n\nReturns the number of all distinct integer partitions of n. The function uses Euler pentagonal number theorem for recursive formula. For more details see OEIS sequence A000041. Note that _numpart(0) = 1 by convention.\n\n\n\n"
},

{
    "location": "ytabs/#Base.conj-Tuple{AbstractAlgebra.Generic.Partition}",
    "page": "Partitions and Young tableaux",
    "title": "Base.conj",
    "category": "method",
    "text": "conj(part::Partition)\n\nReturns the conjugated partition of part, i.e. the partition corresponding to the Young diagram of part reflected through the main diagonal.\n\nExamples:\n\njulia> p = Partition([4,2,1,1,1])\n4₁2₁1₃\n\njulia> conj(p)\n5₁2₁1₂\n\n\n\n"
},

{
    "location": "ytabs/#Base.conj-Tuple{AbstractAlgebra.Generic.Partition,Array{T,1} where T}",
    "page": "Partitions and Young tableaux",
    "title": "Base.conj",
    "category": "method",
    "text": "conj(part::Partition, v::Vector)\n\nReturns the conjugated partition of part together with permuted vector v.\n\n\n\n"
},

{
    "location": "ytabs/#Array-interface-1",
    "page": "Partitions and Young tableaux",
    "title": "Array interface",
    "category": "section",
    "text": "Partition is a concrete subtype of AbstractVector{Int} and implements the following standard Array interface:size(::Generic.Partition)\ngetindex(::Generic.Partition, i::Integer)\nsetindex!(::Generic.Partition, v::Integer, i::Integer)These functions work on the level of p.part vector. Additionally setindex! will try to prevent uses which result in non-valid (i.e. non-decreasing) partition vectors.One can easily iterate over all partitions of n using the AllParts type:Generic.AllPartsThe number all all partitions can be computed by the hidden function _numpart. Much faster implementation is available in Nemo.jl.Generic._numpartSince Partition is a subtype of AbstractVector generic functions which operate on vectors should work in general. However the meaning of conj has been changed to agree with the traditional understanding of conjugation of Partitions:conj(::Generic.Partition)\nconj(::Generic.Partition, v::Vector)"
},

{
    "location": "ytabs/#Young-Diagrams-and-Young-Tableaux-1",
    "page": "Partitions and Young tableaux",
    "title": "Young Diagrams and Young Tableaux",
    "category": "section",
    "text": "Mathematicaly speaking Young diagram is a diagram which consists of rows of square boxes such that the number of boxes in each row is no less than the number of boxes in the previous row. For example partition 4_1 3_2 1 represents the following diagram.┌───┬───┬───┬───┐\n│   │   │   │   │\n├───┼───┼───┼───┘\n│   │   │   │\n├───┼───┼───┤\n│   │   │   │\n├───┼───┴───┘\n│   │\n└───┘Young Tableau is formally a bijection between the set of boxes of a Young Diagram and the set 1 ldots n. If a bijection is increasing along rows and columns of the diagram it is referred to as standard. For example┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┼───┤\n│ 8 │ 9 │10 │\n├───┼───┴───┘\n│11 │\n└───┘is a standard Young tableau of 4_1 3_2 1 where the bijection assigns consecutive natural numbers to consecutive (row-major) cells."
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.YoungTableau",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.YoungTableau",
    "category": "type",
    "text": "YoungTableau(part::Partition[, fill::Vector{Int}=collect(1:sum(part))])  <: AbstractArray{Int, 2}\n\nReturns the Young tableaux of partition part, filled linearly by fill vector. Note that fill vector is in row-major format.Fields:part - the partition defining Young diagram\nfill - the row-major fill vector: the entries of the diagram.\n\nExamples:\n\njulia> p = Partition([4,3,1]); y = YoungTableau(p)\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> y.part\n4₁3₁1₁\n\njulia> y.fill\n8-element Array{Int64,1}:\n 1\n 2\n 3\n 4\n 5\n 6\n 7\n 8\n\n\n\n"
},

{
    "location": "ytabs/#Constructors-1",
    "page": "Partitions and Young tableaux",
    "title": "Constructors",
    "category": "section",
    "text": "In AbstractAlgebra.jl Young tableau are implemented as essentially row-major sparse matrices, i.e. YoungTableau <: AbstractArray{Int,2} but only the defining Partition and the (row-major) fill-vector is stored.Generic.YoungTableauFor convenience there exists an alternative constructor of YoungTableau, which accepts a vector of integers and constructs Partition internally.YoungTableau(p::Vector{Integer}[, fill=collect(1:sum(p))])"
},

{
    "location": "ytabs/#Base.size-Tuple{AbstractAlgebra.Generic.YoungTableau}",
    "page": "Partitions and Young tableaux",
    "title": "Base.size",
    "category": "method",
    "text": "size(Y::YoungTableau)\n\nReturn size of the smallest array containing Y, i.e. the tuple of the number of rows and the number of columns of Y.\n\nExamples:\n\njulia> y = YoungTableau([4,3,1]); size(y)\n(3, 4)\n\n\n\n"
},

{
    "location": "ytabs/#Base.getindex-Tuple{AbstractAlgebra.Generic.YoungTableau,Integer}",
    "page": "Partitions and Young tableaux",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(Y::YoungTableau, n::Integer)\n\nReturn the column-major linear index into the size(Y)-array. If a box is outside of the array return 0.\n\nExamples:\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> y[1]\n1\n\njulia> y[2]\n5\n\njulia> y[4]\n2\n\njulia> y[6]\n0\n\n\n\n"
},

{
    "location": "ytabs/#Base.conj-Tuple{AbstractAlgebra.Generic.YoungTableau}",
    "page": "Partitions and Young tableaux",
    "title": "Base.conj",
    "category": "method",
    "text": "conj(Y::YoungTableau)\n\nReturns the conjugated tableau, i.e. the tableau reflected through the main diagonal.\n\nExamples\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> conj(y)\n┌───┬───┬───┐\n│ 1 │ 5 │ 8 │\n├───┼───┼───┘\n│ 2 │ 6 │\n├───┼───┤\n│ 3 │ 7 │\n├───┼───┘\n│ 4 │\n└───┘\n\n\n\n"
},

{
    "location": "ytabs/#Array-interface-2",
    "page": "Partitions and Young tableaux",
    "title": "Array interface",
    "category": "section",
    "text": "To make YoungTableaux array-like we implement the following functions:size(::Generic.YoungTableau)\ngetindex(::Generic.YoungTableau, n::Integer)Also the double-indexing corresponds to (row, column) access to an abstract array.julia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> y[1,2]\n2\n\njulia> y[2,3]\n7\n\njulia> y[3,2]\n0Functions defined for AbstractArray type based on those (e.g. length) should work. Again, as in the case of Partition the meaning of conj is altered to reflect the usual meaning for Young tableaux:conj(::Generic.YoungTableau)"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.setyoungtabstyle",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.setyoungtabstyle",
    "category": "function",
    "text": "setyoungtabstyle(format::Symbol)\n\nSelect the style in which Young tableaux are displayed (in REPL or in general as string). This can be either:array - as matrices of integers, or\n:diagram - as filled Young diagrams (the default).The difference is purely esthetical.\n\nExamples:\n\njulia> Generic.setyoungtabstyle(:array)\n:array\n\njulia> p = Partition([4,3,1]); YoungTableau(p)\n 1  2  3  4\n 5  6  7\n 8\n\njulia> Generic.setyoungtabstyle(:diagram)\n:diagram\n\njulia> YoungTableau(p)\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\n\n\n"
},

{
    "location": "ytabs/#Pretty-printing-1",
    "page": "Partitions and Young tableaux",
    "title": "Pretty-printing",
    "category": "section",
    "text": "Similarly to permutations we have two methods of displaying Young Diagrams:Generic.setyoungtabstyle"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.matrix_repr-Tuple{AbstractAlgebra.Generic.YoungTableau}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.matrix_repr",
    "category": "method",
    "text": "matrix_repr(Y::YoungTableau)\n\nConstruct sparse integer matrix representing the tableau.\n\nExamples:\n\njulia> y = YoungTableau([4,3,1]);\n\njulia> matrix_repr(y)\n3×4 SparseMatrixCSC{Int64,Int64} with 8 stored entries:\n  [1, 1]  =  1\n  [2, 1]  =  5\n  [3, 1]  =  8\n  [1, 2]  =  2\n  [2, 2]  =  6\n  [1, 3]  =  3\n  [2, 3]  =  7\n  [1, 4]  =  4\n\n\n\n"
},

{
    "location": "ytabs/#Base.fill!-Tuple{AbstractAlgebra.Generic.YoungTableau,AbstractArray{#s1,1} where #s1<:Integer}",
    "page": "Partitions and Young tableaux",
    "title": "Base.fill!",
    "category": "method",
    "text": "fill!(Y::YoungTableaux, V::Vector{<:Integer})\n\nReplace the fill vector Y.fill by V. No check if the resulting tableau is standard (i.e. increasing along rows and columns) is performed.\n\nExamples:\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> fill!(y, [2:9...])\n┌───┬───┬───┬───┐\n│ 2 │ 3 │ 4 │ 5 │\n├───┼───┼───┼───┘\n│ 6 │ 7 │ 8 │\n├───┼───┴───┘\n│ 9 │\n└───┘\n\n\n\n"
},

{
    "location": "ytabs/#Ulitility-functions-1",
    "page": "Partitions and Young tableaux",
    "title": "Ulitility functions",
    "category": "section",
    "text": "matrix_repr(::Generic.YoungTableau)\nfill!(::Generic.YoungTableau, ::AbstractVector{<:Integer})"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.character-Tuple{AbstractAlgebra.Generic.Partition}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.character",
    "category": "method",
    "text": "character(lambda::Partition)\n\nReturn the lambda-th irreducible character of permutation group on sum(lambda) symbols. The returned character function is of the following signature:chi(p::perm[, check::Bool=true]) -> BigIntThe function checks (if p belongs to the appropriate group) can be switched off by calling chi(p, false). The values computed by chi are cached in look-up table.The computation follows the Murnaghan-Nakayama formula: chi_lambda(sigma) = sum_textrimhook xisubset lambda(-1)^ll(lambdabackslashxi) chi_lambda backslashxi(tildesigma) where lambdabackslashxi denotes the skew diagram of lambda with xi removed, ll denotes the leg-length (i.e. number of rows - 1) and tildesigma is permutation obtained from sigma by the removal of the longest cycle.For more details see e.g. Chapter 2.8 of Group Theory and Physics by S.Sternberg.\n\nExamples\n\njulia> G = PermutationGroup(4)\nPermutation group over 4 elements\n\njulia> chi = character(Partition([3,1])) # character of the regular representation\n(::char) (generic function with 2 methods)\n\njulia> chi(G())\n3\n\njulia> chi(perm\"(1,3)(2,4)\")\n-1\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.character-Tuple{AbstractAlgebra.Generic.Partition,perm}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.character",
    "category": "method",
    "text": "character(lambda::Partition, p::perm, check::Bool=true) -> BigInt\n\nReturns the value of lambda-th irreducible character of the permutation group on permutation p.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.character-Tuple{AbstractAlgebra.Generic.Partition,AbstractAlgebra.Generic.Partition}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.character",
    "category": "method",
    "text": "character(lambda::Partition, mu::Partition, check::Bool=true) -> BigInt\n\nReturns the value of lambda-th irreducible character on the conjugacy class represented by partition mu.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.rowlength",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.rowlength",
    "category": "function",
    "text": "rowlength(Y::YoungTableau, i, j)\n\nReturn the row length of Y at box (i,j), i.e. the number of boxes in the i-th row of the diagram of Y located to the right of the (i,j)-th box.\n\nExamples\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> Generic.rowlength(y, 1,2)\n2\n\njulia> Generic.rowlength(y, 2,3)\n0\n\njulia> Generic.rowlength(y, 3,3)\n0\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.collength",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.collength",
    "category": "function",
    "text": "collength(Y::YoungTableau, i, j)\n\nReturn the column length of Y at box (i,j), i.e. the number of boxes in the j-th column of the diagram of Y located below of the (i,j)-th box.\n\nExamples\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> Generic.collength(y, 1,1)\n2\n\njulia> Generic.collength(y, 1,3)\n1\n\njulia> Generic.collength(y, 2,4)\n0\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.hooklength",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.hooklength",
    "category": "function",
    "text": "hooklength(Y::YoungTableau, i, j)\n\nReturn the hook-length of an element in Y at position (i,j), i.e the number of cells in the i-th row to the rigth of (i,j)-th box, plus the number of cells in the j-th column below the (i,j)-th box, plus 1.Return 0 for (i,j) not in the tableau Y.\n\nExamples\n\njulia> y = YoungTableau([4,3,1])\n┌───┬───┬───┬───┐\n│ 1 │ 2 │ 3 │ 4 │\n├───┼───┼───┼───┘\n│ 5 │ 6 │ 7 │\n├───┼───┴───┘\n│ 8 │\n└───┘\n\njulia> hooklength(y, 1,1)\n6\n\njulia> hooklength(y, 1,3)\n3\n\njulia> hooklength(y, 2,4)\n0\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.dim-Tuple{AbstractAlgebra.Generic.YoungTableau}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.dim",
    "category": "method",
    "text": "dim(Y::YoungTableau) -> BigInt\n\nReturns the dimension (using hook-length formula) of the irreducible representation of permutation group S_n associated the partition Y.part.Since the computation overflows easily BigInt is returned. You may perform the computation of the dimension in different type by calling dim(Int, Y).\n\nExamples\n\njulia> dim(YoungTableau([4,3,1]))\n70\n\njulia> dim(YoungTableau([3,1])) # the regular representation of S_4\n3\n\n\n\n"
},

{
    "location": "ytabs/#Characters-of-permutation-grups-1",
    "page": "Partitions and Young tableaux",
    "title": "Characters of permutation grups",
    "category": "section",
    "text": "Irreducible characters (at least over field of characteristic 0) of the full group of permutations S_n correspond via Specht modules to partitions of n.character(::Generic.Partition)\ncharacter(lambda::Generic.Partition, p::Generic.perm)\ncharacter(lambda::Generic.Partition, mu::Generic.Partition)The values computed by characters are cached in an internal dictionary Dict{Tuple{BitVector,Vector{Int}}, BigInt}. Note that all of the above functions return BigInts. If you are sure that the computations do not overflow, variants of the last two functions using Int are available:character(::Type{Int}, lambda::Partition, p::perm[, check::Bool=true])\ncharacter(::Type{Int}, lambda::Partition, mu::Partition[, check::Bool=true])The dimension dim lambda of the irreducible module corresponding to partition lambda can be computed using Hook length formulaGeneric.rowlength\nGeneric.collength\nhooklength\ndim(::Generic.YoungTableau)The the character associated with Y.part can also be used to compute the dimension, but as it is expected the Murnaghan-Nakayama is much slower even though (due to caching) consecutive calls are fast:julia> λ = Partition(collect(12:-1:1))\n12₁11₁10₁9₁8₁7₁6₁5₁4₁3₁2₁1₁\n\njulia> @time dim(YoungTableau(λ))\n  0.224430 seconds (155.77 k allocations: 7.990 MiB)\n9079590132732747656880081324531330222983622187548672000\n\njulia> @time dim(YoungTableau(λ))\n  0.000038 seconds (335 allocations: 10.734 KiB)\n9079590132732747656880081324531330222983622187548672000\n\njulia> G = PermutationGroup(sum(λ))\nPermutation group over 78 elements\n\njulia> @time character(λ, G())\n 24.154105 seconds (58.13 M allocations: 3.909 GiB, 42.84% gc time)\n9079590132732747656880081324531330222983622187548672000\n\njulia> @time character(λ, G())\n  0.001439 seconds (195 allocations: 24.453 KiB)\n9079590132732747656880081324531330222983622187548672000"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.partitionseq",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.partitionseq",
    "category": "function",
    "text": "partitionseq(lambda::Partition)\n\nReturns a sequence (as BitVector) of falses and trues constructed from lambda: tracing the lower contour of the Young Diagram associated to lambda from left to right a true is inserted for every horizontal and false for every vertical step. The sequence always starts with true and ends with false.\n\n\n\npartitionseq(seq::BitVector)\n\nReturns the essential part of the sequence seq, i.e. a subsequence starting at first true and ending at last false.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.isrimhook-Tuple{BitArray{1},Int64,Int64}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.isrimhook",
    "category": "method",
    "text": "isrimhook(R::BitVector, idx::Int, len::Int)\n\nR[idx:idx+len] forms a rim hook in the Young Diagram of partition corresponding to R iff R[idx] == true and R[idx+len] == false.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.MN1inner",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.MN1inner",
    "category": "function",
    "text": "MN1inner(R::BitVector, mu::Partition, t::Integer, charvals)\n\nReturns the value of lambda-th irreducible character on conjugacy class of permutations represented by partition mu, where R is the (binary) partition sequence representing lambda. Values already computed are stored in charvals::Dict{Tuple{BitVector,Vector{Int}}, Int}. This is an implementation (with slight modifications) of the Murnaghan-Nakayama formula as described inDan Bernstein,\n\"The computational complexity of rules for the character table of Sn\"\n_Journal of Symbolic Computation_, 37(6), 2004, p. 727-748.\n\n\n\n"
},

{
    "location": "ytabs/#Low-level-functions-and-characters-1",
    "page": "Partitions and Young tableaux",
    "title": "Low-level functions and characters",
    "category": "section",
    "text": "As mentioned above character functions use the Murnaghan-Nakayama rule for evaluation. The implementation followsDan Bernstein, The computational complexity of rules for the character table of S_n Journal of Symbolic Computation, 37 (6), 2004, p. 727-748,implementing the following functions. For precise definitions and meaning please consult the paper cited.Generic.partitionseq\nisrimhook(::BitVector, ::Int, ::Int)\nGeneric.MN1inner"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.SkewDiagram",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.SkewDiagram",
    "category": "type",
    "text": "SkewDiagram(lambda::Partition, mu::Partition) <: AbstractArray{Int, 2}\n\nImplements a skew diagram, i.e. a difference of two Young diagrams represented by partitions lambda and mu. (below dots symbolise the removed entries)\n\nExamples\n\njulia> l = Partition([4,3,2])\n4₁3₁2₁\n\njulia> m = Partition([3,1,1])\n3₁1₂\n\njulia> xi = SkewDiagram(l,m)\n3×4 AbstractAlgebra.Generic.SkewDiagram:\n ⋅  ⋅  ⋅  1\n ⋅  1  1\n ⋅  1\n\n\n\n\n"
},

{
    "location": "ytabs/#Base.size-Tuple{AbstractAlgebra.Generic.SkewDiagram}",
    "page": "Partitions and Young tableaux",
    "title": "Base.size",
    "category": "method",
    "text": "size(xi::SkewDiagram)\n\nReturn the size of array where xi is minimally contained. See size(Y::YoungTableau) for more details.\n\n\n\n"
},

{
    "location": "ytabs/#Base.in-Union{Tuple{T}, Tuple{Tuple{T,T},SkewDiagram}} where T<:Integer",
    "page": "Partitions and Young tableaux",
    "title": "Base.in",
    "category": "method",
    "text": "in(t::Tuple{T,T}, xi::SkewDiagram) where T<:Integer\n\nChecks if box at position (i,j) belongs to the skew diagram xi.\n\n\n\n"
},

{
    "location": "ytabs/#Base.getindex-Tuple{AbstractAlgebra.Generic.SkewDiagram,Integer}",
    "page": "Partitions and Young tableaux",
    "title": "Base.getindex",
    "category": "method",
    "text": "getindex(xi::SkewDiagram, n::Integer)\n\nReturn 1 if linear index n corresponds to (column-major) entry in xi.lam which is not contained in xi.mu. Otherwise return 0.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.isrimhook-Tuple{AbstractAlgebra.Generic.SkewDiagram}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.isrimhook",
    "category": "method",
    "text": "isrimhook(xi::SkewDiagram)\n\nChecks if xi represents a rim-hook diagram, i.e. its diagram is edge-connected and contains no 2times 2 squares.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.leglength",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.leglength",
    "category": "function",
    "text": "leglength(xi::SkewDiagram[, check::Bool=true])\n\nComputes the leglength of a rim-hook xi, i.e. the number of rows with non-zero entries minus one. If check is false function will not check whether xi is actually a rim-hook.\n\n\n\n"
},

{
    "location": "ytabs/#AbstractAlgebra.Generic.matrix_repr-Tuple{AbstractAlgebra.Generic.SkewDiagram}",
    "page": "Partitions and Young tableaux",
    "title": "AbstractAlgebra.Generic.matrix_repr",
    "category": "method",
    "text": "matrix_repr(xi::SkewDiagram)\n\nReturns a sparse representation of the diagram xi, i.e. a sparse array A where A[i,j] == 1 if and only if (i,j) is in xi.lam but not in xi.mu.\n\n\n\n"
},

{
    "location": "ytabs/#Skew-Diagrams-1",
    "page": "Partitions and Young tableaux",
    "title": "Skew Diagrams",
    "category": "section",
    "text": "Skew diagrams are formally differences of two Young diagrams. Given lambda and mu, two partitions of n+m and m (respectively). Suppose that each of cells of mu is a cell of lambda (i.e. parts of mu are no greater than the corresponding parts of lambda). Then the skew diagram denoted by lambdamu is the set theoretic difference the of sets of boxes, i.e. is a diagram with exactly n boxes:Generic.SkewDiagramSkewDiagram implements array interface with the following functions:size(xi::Generic.SkewDiagram)\nin(t::Tuple{T,T}, xi::Generic.SkewDiagram) where T<:Integer\ngetindex(xi::Generic.SkewDiagram, n::Integer)The support for skew diagrams is very rudimentary. The following functions are available:isrimhook(::Generic.SkewDiagram)\nleglength\nmatrix_repr(::Generic.SkewDiagram)"
},

{
    "location": "matrix_spaces/#",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "page",
    "text": ""
},

{
    "location": "matrix_spaces/#Matrix-Interface-1",
    "page": "Matrix Interface",
    "title": "Matrix Interface",
    "category": "section",
    "text": "Generic matrices are supported in AbstractAlgebra.jl. Both the space of mtimes n  matrices and the algebra (ring) of mtimes m matrices are supported.As the space of mtimes n matrices over a commutative ring is not itself a commutative ring, not all of the Ring interface needs to be implemented for such matrices in.In particular, the following functions do not need to be implemented: isdomain_type, needs_parentheses, displayed_with_minus_in_front, show_minus_one and divexact. The canonical_unit function should be implemented, but simply needs to return the corresponding value for entry 1 1 (the function is never called on empty matrices).For matrix algebras, all of the ring interface must be implemented.Note that AbstractAlgebra.jl matrices are not the same as Julia matrices. We store a base ring in our matrix and matrices are row major instead of column major in order to support the numerous large C libraries that use this convention.All AbstractAlgebra.jl matrices are assumed to be mutable. This is usually critical to performance."
},

{
    "location": "matrix_spaces/#Types-and-parents-1",
    "page": "Matrix Interface",
    "title": "Types and parents",
    "category": "section",
    "text": "AbstractAlgebra provides two abstract types for matrix spaces and their elements:MatSpace{T} is the abstract type for matrix space parent types\nMatElem{T} is the abstract type for matrix types belonging to a matrix spaceIt also provides two abstract types for matrix algebras and their elements:MatAlgebra{T} is the abstract type for matrix algebra parent types\nMatAlgElem{T} is the abstract type for matrix types belonging to a matrix algebraNote that these abstract types are parameterised. The type T should usually be the type of elements of the matrices.Matrix spaces and matrix algebras should be made unique on the system by caching parent objects (unless an optional cache parameter is set to false). Matrix spaces and algebras should at least be distinguished based on their base (coefficient) ring and the  dimensions of the matrices in the space.See src/generic/GenericTypes.jl for an example of how to implement such a cache (which usually makes use of a dictionary)."
},

{
    "location": "matrix_spaces/#Required-functionality-for-matrices-1",
    "page": "Matrix Interface",
    "title": "Required functionality for matrices",
    "category": "section",
    "text": "In addition to the required (relevant) functionality for the Ring interface (see above), the following functionality is required for the Matrix interface.We suppose that R is a fictitious base ring (coefficient ring) and that S is a space of mtimes n matrices over R, or algebra of mtimes m matrices with parent object S of type MyMatSpace{T} or MyMatAlgebra{T}, respectively. We also assume the matrices in the space have type MyMat{T}, where T is the type of elements of the base (element) ring.Of course, in practice these types may not be parameterised, but we use parameterised types here to make the interface clearer.Note that the type T must (transitively) belong to the abstract type RingElem.Currently only matrices over commutative rings are supported."
},

{
    "location": "matrix_spaces/#Constructors-1",
    "page": "Matrix Interface",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors, the following constructors, taking an array of elements, must be available.(S::MyMatSpace{T})(A::Array{T, 2}) where T <: AbstractAlgebra.RingElem\n(S::MyMatAlgebra{T})(A::Array{T, 2}) where T <: AbstractAlgebra.RingElemCreate the matrix in the given space/algebra whose (i j) entry is given by A[i, j].(S::MyMatSpace{T})(A::Array{S, 2}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}\n(S::MyMatAlgebra{T})(A::Array{S, 2}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}Create the matrix in the given space/algebra whose (i j) entry is given by A[i, j], where S is the type of elements that can be coerced into the base ring of the matrix.(S::MyMatSpace{T})(A::Array{S, 1}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}\n(S::MyMatAlgebra{T})(A::Array{S, 1}) where {S <: AbstractAlgebra.RingElem, T <: AbstractAlgebra.RingElem}Create the matrix in the given space/algebra of matrices (with dimensions mtimes n say), whose (i j) entry is given by A[i*(n - 1) + j] and where S is the type of elements that can be coerced into the base ring of the matrix.ExamplesS = MatrixSpace(QQ, 2, 3)\nT = MatrixAlgebra(QQ, 2)\n\nM1 = S(Rational{BigInt}[2 3 1; 1 0 4])\nM2 = S(BigInt[2 3 1; 1 0 4])\nM3 = S(BigInt[2, 3, 1, 1, 0, 4])\n\nN1 = T(Rational{BigInt}[2 3; 1 0])\nN2 = T(BigInt[2 3; 1 0])\nN3 = T(BigInt[2, 3, 1, 1])It is also possible to create matrices (in a matrix space only) directly, without first creating the corresponding matrix space (the inner constructor being called directly). Note that to support this, matrix space parent objects don\'t contain a reference to their parent. Instead, parents are constructed on-the-fly if requested. (The same strategy is used for matrix algebras.)matrix(R::Ring, arr::Array{T, 2}) where T <: AbstractAlgebra.RingElemGiven an mtimes n Julia matrix of entries, construct the corresponding AbstractAlgebra.jl matrix over the given ring R, assuming all the entries can be coerced into R.matrix(R::Ring, r::Int, c::Int, A::Array{T, 1}) where T <: AbstractAlgebra.RingElemConstruct the given rtimes c AbstractAlgebra.jl matrix over the ring R whose (i j) entry is given by A[c*(i - 1) + j], assuming that all the entries can be coerced into R.zero_matrix(R::Ring, r::Int, c::Int)Construct the rtimes c AbstractAlgebra.jl zero matrix over the ring R.identity_matrix(R::Ring, n::Int)Construct the ntimes n AbstractAlgebra.jl identity matrix over the ring R.The following functions are available for matrices in both matrix algebras and matrix spaces.similar(x::MyMat{T}) where T <: AbstractAlgebra.RingElemConstruct the zero matrix with the same dimensions and base ring as the given matrix.similar(x::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElemConstruct the rtimes c zero matrix with the same base ring as the given matrix. If x belongs to a matrix algebra and r neq c, an exception is raised.ExamplesM = matrix(ZZ, BigInt[3 1 2; 2 0 1])\nN = matrix(ZZ, 3, 2, BigInt[3, 1, 2, 2, 0, 1])\nP = zero_matrix(ZZ, 3, 2)\nQ = identity_matrix(ZZ, 4)\nC = similar(P)\nD = similar(Q, 4, 5)\n\nR = MatrixAlgebra(ZZ, 2)\nM = R()\nF = similar(M)"
},

{
    "location": "matrix_spaces/#Basic-manipulation-of-matrices-1",
    "page": "Matrix Interface",
    "title": "Basic manipulation of matrices",
    "category": "section",
    "text": "rows(f::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the number of rows of the given matrix.cols(f::MyMat{T}) where T <: AbstractAlgebra.RingElemReturns the number of columns of the given matrix.getindex(M::MyMat{T}, r::Int, c::Int) where T <: AbstractAlgebra.RingElemReturn the (i j)-th entry of the matrix M.setindex!(M::MyMat{T}, d::T, r::Int, c::Int) where T <: AbstractAlgebra.RingElemSet the (i j)-th entry of the matrix M to d, which is assumed to be in the base ring of the matrix. The matrix must have such an entry and the matrix is mutated in place and not returned from the function.ExamplesM = matrix(ZZ, BigInt[2 3 0; 1 1 1])\n\nm = rows(M)\nn = cols(M)\nM[1, 2] = BigInt(4)\nc = M[1, 1]"
},

{
    "location": "matrix_spaces/#Transpose-1",
    "page": "Matrix Interface",
    "title": "Transpose",
    "category": "section",
    "text": "transpose(::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the transpose of the given matrix.The standard Julia tick notation can also be used for transposing a matrix.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nB = transpose(A)\nC = A\'"
},

{
    "location": "matrix_spaces/#Optional-functionality-for-matrices-1",
    "page": "Matrix Interface",
    "title": "Optional functionality for matrices",
    "category": "section",
    "text": "Especially when wrapping C libraries, some functions are best implemented directly, rather than relying on the generic functionality. The following are all provided by the AbstractAlgebra.jl generic code, but can optionally be implemented directly for performance reasons."
},

{
    "location": "matrix_spaces/#Optional-constructors-1",
    "page": "Matrix Interface",
    "title": "Optional constructors",
    "category": "section",
    "text": "The following can only be used to construct matrices in a matrix space, not a matrix algebra.eye(M::MyMat{T}) where T <: AbstractAlgebra.RingElemConstruct the identity matrix with the same dimensions and base ring as the given matrix.eye(M::MyMat{T}, n::Int) where T <: AbstractAlgebra.RingElemConstruct the ntimes n identity matrix with the same base ring as the given matrix.ExamplesM = matrix(ZZ, BigInt[1 2 3; 4 5 6])\n\nN = eye(M)\nP = eye(M, 2)"
},

{
    "location": "matrix_spaces/#Optional-submatrices-1",
    "page": "Matrix Interface",
    "title": "Optional submatrices",
    "category": "section",
    "text": "The following are only available for matrix spaces, not for matrix algebras.sub(M::MyMat{T}, rows::UnitRange{Int}, cols::UnitRange{Int}) where T <: AbstractAlgebra.RingElemReturn a new matrix with the same entries as the submatrix with the given range of rows and columns.ExamplesM = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])\n\nN1 = M[1:2, :]\nN2 = M[:, :]\nN3 = M[2:3, 2:3]"
},

{
    "location": "matrix_spaces/#Optional-row-swapping-1",
    "page": "Matrix Interface",
    "title": "Optional row swapping",
    "category": "section",
    "text": "swap_rows!(M::MyMat{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElemSwap the rows of M in place. The function does not return the mutated matrix (since matrices are assumed to be mutable in AbstractAlgebra.jl).ExamplesM = identity_matrix(ZZ, 3)\n\nswap_rows!(M, 1, 2)"
},

{
    "location": "matrix_spaces/#Optional-concatenation-1",
    "page": "Matrix Interface",
    "title": "Optional concatenation",
    "category": "section",
    "text": "The following are only available for matrix spaces, not for matrix algebras.hcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the horizontal concatenation of M and N. It is assumed that the number of rows of M and N are the same.vcat(M::MyMat{T}, N::MyMat{T}) where T <: AbstractAlgebra.RingElemReturn the vertical concatenation of M and N. It is assumed that the number of columns of M and N are the same.ExamplesM = matrix(ZZ, BigInt[1 2 3; 2 3 4; 3 4 5])\nN = matrix(ZZ, BigInt[1 0 1; 0 1 0; 1 0 1])\n\nP = hcat(M, N)\nQ = vcat(M, N)"
},

{
    "location": "matrix/#",
    "page": "Generic matrices",
    "title": "Generic matrices",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "matrix/#Generic-matrices-1",
    "page": "Generic matrices",
    "title": "Generic matrices",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of dense matrices over any computable commutative ring R. Generic matrices over a commutative ring are implemented in src/generic/Matrix.jl. Much of the functionality there covers both matrix spaces and matrix algebras.Functions specific to generic matrix algebras of mtimes m matrices are implemented in  src/generic/MatrixAlgebra.jl.As well as implementing the entire Matrix interface, including the optional functionality, there are many additional generic algorithms implemented for matrix spaces. We describe this functionality below.All of this generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names of functions."
},

{
    "location": "matrix/#Types-and-parent-objects-1",
    "page": "Generic matrices",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Generic matrices in AbstractAlgebra.jl have type Generic.Mat{T} for matrices in a matrix space, or Generic.MatAlgElem{T} for matrices in a matrix algebra, where T is the type of elements of the matrix. Internally, generic matrices are implemented using an object wrapping a Julia two dimensional array, though they are not themselves Julia arrays. See the file src/generic/GenericTypes.jl for details.Parents of generic matrices (matrix spaces) have type Generic.MatSpace{T}. Parents of matrices in a matrix algebra have type Generic.MatAlgebra{T}.The generic matrix types (matrix spaces) belong to the abstract type AbstractAlgebra.MatElem{T} and the matrix space parent types belong to AbstractAlgebra.MatSpace{T}. Similarly the generic matrix algebra matrix types belong to the abstract type AbstractAlgebra.MatAlgElem{T} and the parent types belong to  AbstractAlgebra.MatAlgebra{T} Note that both the concrete type of a matrix space parent object and the abstract class it belongs to have the name MatElem, therefore disambiguation is required to specify which is intended. The same is true for the abstract types for matrix spaces and their elements.The dimensions and base ring R of a generic matrix are stored in its parent object, however to allow creation of matrices without first creating the matrix space parent, generic matrices in Julia do not contain a reference to their parent. They contain the row and column numbers (or degree, in the case of matrix algebras) and the base ring on a per matrix basis. The parent object can then be reconstructed from this data on demand."
},

{
    "location": "matrix/#Matrix-space-constructors-1",
    "page": "Generic matrices",
    "title": "Matrix space constructors",
    "category": "section",
    "text": "A matrix space in AbstractAlgebra.jl represents a collection of all matrices with given dimensions and base ring.In order to construct matrices in AbstractAlgebra.jl, one can first constructs the matrix space itself. This is accomplished with the following constructor. We discuss creation of matrix algebras separately in a dedicated section elsewhere in the documentation.MatrixSpace(R::Ring, rows::Int, cols::Int; cache::Bool=true)Construct the space of matrices with the given number of rows and columns over the given base ring. By default such matrix spaces are cached based on the base ring and numbers of rows and columns. If the optional named parameter cached is set to false, no caching occurs.Here are some examples of creating matrix spaces and making use of the resulting parent objects to coerce various elements into the matrix space.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S()\nB = S(12)\nC = S(R(11))We also allow matrices over a given base ring to be constructed directly (see the Matrix interface)."
},

{
    "location": "matrix/#Matrix-element-constructors-1",
    "page": "Generic matrices",
    "title": "Matrix element constructors",
    "category": "section",
    "text": "In addition to coercing elements into a matrix space as above, we provide the following functions for constructing explicit matrices.Also see the Matrix interface for a list of other ways to create matrices.R[a b c...;...]Create the matrix over the base ring R consisting of the given rows (separated by semicolons). Each entry is coerced into R  automatically. Note that parentheses may be placed around individual entries if the lists would otherwise be ambiguous, e.g.  R[1 2; 2 (-3)].Beware that this syntax does not support the creation of column vectors. See the notation below for creating those.R[a b c...]Create the row vector with entries in R consisting of the given entries (separated by spaces). Each entry is coerced into R automatically. Note that parentheses may be placed around individual entries if the list would otherwise be ambiguous, e.g. R[1 2 (-3)].R[a b c...]\'Create the column vector with entries in R consisting of the given entries (separated by spaces). Each entry is coerced into R automatically. Observe the dash that is used to transpose the row vector notation (for free) to turn it into a column vector. Note that parentheses may be placed around individual entries if the list would otherwise be ambiguous, e.g. R[1 2 (-3)]\'.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nM = R[t + 1 1; t^2 0]\nN = R[t + 1 2 t]\nP = R[1 2 t]\'"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.sub-Tuple{MatElem,Int64,Int64,Int64,Int64}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.sub",
    "category": "method",
    "text": "sub(M::AbstractAlgebra.MatElem, r1::Int, c1::Int, r2::Int, c2::Int)\n\nReturn a copy of the submatrix of M from (r1 c1) to (r2 c2) inclusive. Note that is the copy is modified, the original matrix is not.\n\n\n\n"
},

{
    "location": "matrix/#Submatrices-1",
    "page": "Generic matrices",
    "title": "Submatrices",
    "category": "section",
    "text": "In addition to the functionality described in the Matrix interface for taking submatrices of a matrix, the following function variant is also available.sub(::MatElem, ::Int, ::Int, ::Int, ::Int)ExamplesM = ZZ[1 2 3; 2 3 4]\n\nN = sub(M, 1, 1, 2, 2)"
},

{
    "location": "matrix/#Matrix-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic matrices",
    "title": "Matrix functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "Most of the following generic functionality is available for both matrix spaces and matrix algebras. Exceptions include functions that do not return or accept square matrices or which cannot specify a parent. Such functions include solve and nullspace which can\'t be provided for matrix algebras.For details on functionality that is provided for matrix algebras only, see the dedicated section of the documentation."
},

{
    "location": "matrix/#AbstractAlgebra.Generic.rows-Tuple{MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rows",
    "category": "method",
    "text": "rows(a::Generic.MatrixElem)\n\nReturn the number of rows of the given matrix.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.cols-Tuple{MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.cols",
    "category": "method",
    "text": "cols(a::Generic.MatrixElem)\n\nReturn the number of columns of the given matrix.\n\n\n\n"
},

{
    "location": "matrix/#Basic-matrix-functionality-1",
    "page": "Generic matrices",
    "title": "Basic matrix functionality",
    "category": "section",
    "text": "As well as the Ring and Matrix interfaces, the following functions are provided to manipulate matrices and to set and retrieve entries and other basic data associated with the matrices.rows(::MatElem)cols(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\nB = S([R(2) R(3) R(1); t t + 1 t + 2; R(-1) t^2 t^3])\n\nr = rows(B)\nc = cols(B)\nM = A + B\nN = 2 + A\nM1 = deepcopy(A)\nA != B\nisone(one(S)) == true\nV = A[1:2, :]\nW = A^3\nZ = divexact(2*A, 2)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.powers-Tuple{MatElem,Int64}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.powers",
    "category": "method",
    "text": "powers(a::Generic.MatrixElem, d::Int)\n\nReturn an array of matrices M wher Mi + 1 = a^i for i = 0d\n\n\n\n"
},

{
    "location": "matrix/#Powering-1",
    "page": "Generic matrices",
    "title": "Powering",
    "category": "section",
    "text": "powers(::MatElem, ::Int)ExamplesM = ZZ[1 2 3; 2 3 4; 4 5 5]\n\nA = powers(M, 4)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.gram-Tuple{MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.gram",
    "category": "method",
    "text": "gram(x::AbstractAlgebra.MatElem)\n\nReturn the Gram matrix of x, i.e. if x is an rtimes c matrix return the rtimes r matrix whose entries i j are the dot products of the i-th and j-th rows, respectively.\n\n\n\n"
},

{
    "location": "matrix/#Gram-matrix-1",
    "page": "Generic matrices",
    "title": "Gram matrix",
    "category": "section",
    "text": "gram(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nB = gram(A)"
},

{
    "location": "matrix/#LinearAlgebra.tr-Tuple{MatElem}",
    "page": "Generic matrices",
    "title": "LinearAlgebra.tr",
    "category": "method",
    "text": "tr(x::Generic.MatrixElem)\n\nReturn the trace of the matrix a, i.e. the sum of the diagonal elements. We require the matrix to be square.\n\n\n\n"
},

{
    "location": "matrix/#Trace-1",
    "page": "Generic matrices",
    "title": "Trace",
    "category": "section",
    "text": "tr(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nb = tr(A)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.content-Tuple{MatElem}",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.content",
    "category": "method",
    "text": "content(x::Generic.MatrixElem)\n\nReturn the content of the matrix a, i.e. the greatest common divisor of all its entries, assuming it exists.\n\n\n\n"
},

{
    "location": "matrix/#Content-1",
    "page": "Generic matrices",
    "title": "Content",
    "category": "section",
    "text": "content(::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nb = content(A)"
},

{
    "location": "matrix/#Base.:*-Tuple{perm,MatElem}",
    "page": "Generic matrices",
    "title": "Base.:*",
    "category": "method",
    "text": "*(P::Generic.perm, x::Generic.MatrixElem)\n\nApply the pemutation P to the rows of the matrix x and return the result.\n\n\n\n"
},

{
    "location": "matrix/#Permutation-1",
    "page": "Generic matrices",
    "title": "Permutation",
    "category": "section",
    "text": "*(::perm, ::MatElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixSpace(R, 3, 3)\nG = PermGroup(3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\nP = G([1, 3, 2])\n\nB = P*A"
},

{
    "location": "matrix/#LinearAlgebra.lu-Union{Tuple{T}, Tuple{MatElem{T},PermGroup}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.lu",
    "category": "method",
    "text": "lu(A::Generic.MatrixElem{T}, P = PermGroup(rows(A))) where {T <: FieldElement}\n\nReturn a tuple r p L U consisting of the rank of A, a permutation p of A belonging to P, a lower triangular matrix L and an upper triangular matrix U such that p(A) = LU, where p(A) stands for the matrix whose rows are the given permutation p of the rows of A.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.fflu-Union{Tuple{T}, Tuple{MatElem{T},PermGroup}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.fflu",
    "category": "method",
    "text": "fflu(A::Generic.MatrixElem{T}, P = PermGroup(rows(A))) where {T <: RingElement}\n\nReturn a tuple r d p L U consisting of the rank of A, a denominator d, a permutation p of A belonging to P, a lower triangular matrix L and an upper triangular matrix U such that p(A) = LD^1U, where p(A) stands for the matrix whose rows are the given permutation p of the rows of A and such that D is the diagonal matrix diag(p_1 p_1p_2 ldots p_n-2p_n-1 p_n-1 where the p_i are the inverses of the diagonal entries of U. The denominator d is set to pm mboxdet(S) where S is an appropriate submatrix of A (S = A if A is square) and the sign is decided by the parity of the permutation.\n\n\n\n"
},

{
    "location": "matrix/#LU-factorisation-1",
    "page": "Generic matrices",
    "title": "LU factorisation",
    "category": "section",
    "text": "lu{T <: FieldElem}(::MatElem{T}, ::PermGroup)fflu{T <: RingElem}(::MatElem{T}, ::PermGroup)ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 - 2 a - 1 2a])\n\nr, P, L, U = lu(A)\nr, d, P, L, U = fflu(A)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.rref-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rref",
    "category": "method",
    "text": "rref(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturns a tuple (r d A) consisting of the rank r of M and a denominator d in the base ring of M and a matrix A such that Ad is the reduced row echelon form of M. Note that the denominator is not usually minimal.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.rref-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.rref",
    "category": "method",
    "text": "rref(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturns a tuple (r d A) consisting of the rank r of M and a denominator d in the base ring of M and a matrix A such that Ad is the reduced row echelon form of M. Note that the denominator is not usually minimal.\n\n\n\nrref(M::Generic.MatrixElem{T}) where {T <: FieldElement}\n\nReturns a tuple (r A) consisting of the rank r of M and a reduced row echelon form A of M.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.isrref-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.isrref",
    "category": "method",
    "text": "isrref(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.isrref-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.isrref",
    "category": "method",
    "text": "isrref(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\nisrref(M::Generic.MatrixElem{T}) where {T <: FieldElement}\n\nReturn true if M is in reduced row echelon form, otherwise return false.\n\n\n\n"
},

{
    "location": "matrix/#Reduced-row-echelon-form-1",
    "page": "Generic matrices",
    "title": "Reduced row-echelon form",
    "category": "section",
    "text": "rref{T <: RingElem}(::MatElem{T})\nrref{T <: FieldElem}(::MatElem{T})isrref{T <: RingElem}(::MatElem{T})\nisrref{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nM = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n   \nr, A = rref(M)\nisrref(A)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\n\nM = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\n\nr, d, A = rref(M)\nisrref(A)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.hnf-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.hnf",
    "category": "method",
    "text": "hnf(A::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn the upper right row Hermite normal form of A.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.hnf_with_trafo-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.hnf_with_trafo",
    "category": "method",
    "text": "hnf_with_trafo(A)\n\nReturn the tuple H U consisting of the upper right row Hermite normal form H of A together with invertible matrix U such that UA = H.\n\n\n\n"
},

{
    "location": "matrix/#Hermite-normal-form-1",
    "page": "Generic matrices",
    "title": "Hermite normal form",
    "category": "section",
    "text": "hnf{T <: RingElem}(::MatElem{T})\nhnf_with_trafo{T <: RingElem}(::MatElem{T})"
},

{
    "location": "matrix/#LinearAlgebra.det-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.det",
    "category": "method",
    "text": "det(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\n"
},

{
    "location": "matrix/#LinearAlgebra.det-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.det",
    "category": "method",
    "text": "det(M::Generic.MatrixElem{T}) where {T <: FieldElement}\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\ndet(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn the determinant of the matrix M. We assume M is square.\n\n\n\n"
},

{
    "location": "matrix/#Determinant-1",
    "page": "Generic matrices",
    "title": "Determinant",
    "category": "section",
    "text": "det{T <: RingElem}(::MatElem{T})\ndet{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nd = det(A)"
},

{
    "location": "matrix/#LinearAlgebra.rank-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.rank",
    "category": "method",
    "text": "rank(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn the rank of the matrix M.\n\n\n\n"
},

{
    "location": "matrix/#LinearAlgebra.rank-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.rank",
    "category": "method",
    "text": "rank(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturn the rank of the matrix M.\n\n\n\nrank(M::Generic.MatrixElem{T}) where {T <: FieldElement}\n\nReturn the rank of the matrix M.\n\n\n\n"
},

{
    "location": "matrix/#Rank-1",
    "page": "Generic matrices",
    "title": "Rank",
    "category": "section",
    "text": "rank{T <: RingElem}(::MatElem{T})\nrank{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n   \nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nd = rank(A)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.solve-Union{Tuple{T}, Tuple{MatElem{T},MatElem{T}}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve",
    "category": "method",
    "text": "solve(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}\n\nGiven a non-singular ntimes n matrix over a field and an ntimes m matrix over the same field, return x an ntimes m matrix x such that Ax = b. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.solve_rational-Union{Tuple{T}, Tuple{MatElem{T},MatElem{T}}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve_rational",
    "category": "method",
    "text": "solve_rational(M::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}) where T <: RingElement\n\nGiven a non-singular ntimes n matrix over a ring and an ntimes m matrix over the same ring, return a tuple x d consisting of an ntimes m matrix x and a denominator d such that Ax = db. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.solve_triu-Union{Tuple{T}, Tuple{MatElem{T},MatElem{T},Bool}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.solve_triu",
    "category": "method",
    "text": "solve_triu(U::AbstractAlgebra.MatElem{T}, b::AbstractAlgebra.MatElem{T}, unit::Bool = false) where {T <: FieldElement}\n\nGiven a non-singular ntimes n matrix over a field which is upper triangular, and an ntimes m matrix over the same field, return an ntimes m matrix x such that Ax = b. If A is singular an exception is raised. If unit is true then U is assumed to have ones on its diagonal, and the diagonal will not be read.\n\n\n\n"
},

{
    "location": "matrix/#Linear-solving-1",
    "page": "Generic matrices",
    "title": "Linear solving",
    "category": "section",
    "text": "solve{T <: FieldElem}(::MatElem{T}, ::MatElem{T})solve_rational{T <: RingElem}(::MatElem{T}, ::MatElem{T})solve_triu{T <: FieldElem}(::MatElem{T}, ::MatElem{T}, ::Bool)ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\nU = MatrixSpace(K, 3, 1)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\nb = U([2a a + 1 (-a - 1)]\')\n\nx = solve(A, b)\n\nA = S([a + 1 2a + 3 a^2 + 1; K(0) a^2 - 1 2a; K(0) K(0) a])\nb = U([2a a + 1 (-a - 1)]\')\n\nx = solve_triu(A, b, false)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\nU = MatrixSpace(R, 3, 2)\n\nA = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\nb = U([2x x + 1 (-x - 1); x + 1 (-x) x^2]\')\n\nx, d = solve_rational(A, b)\n\nS = MatrixSpace(ZZ, 3, 3)\nT = MatrixSpace(ZZ, 3, 1)\n\nA = S([BigInt(2) 3 5; 1 4 7; 9 2 2])   \nB = T([BigInt(4), 5, 7])\n\nX, d = solve_rational(A, B)"
},

{
    "location": "matrix/#Base.inv-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nGiven a non-singular ntimes n matrix over a ring the tuple X d consisting of an ntimes n matrix X and a denominator d such that AX = dI_n, where I_n is the ntimes n identity matrix. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix/#Base.inv-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "Base.inv",
    "category": "method",
    "text": "inv(M::Generic.MatrixElem{T}) where {T <: RingElement}\n\nGiven a non-singular ntimes n matrix over a ring the tuple X d consisting of an ntimes n matrix X and a denominator d such that AX = dI_n, where I_n is the ntimes n identity matrix. The denominator will be the determinant of A up to sign. If A is singular an exception is raised.\n\n\n\ninv(M::Generic.MatrixElem{T}) where {T <: FieldElement}\n\nGiven a non-singular ntimes n matrix over a field, return an ntimes n matrix X such that AX = I_n where I_n is the ntimes n identity matrix. If A is singular an exception is raised.\n\n\n\n"
},

{
    "location": "matrix/#Inverse-1",
    "page": "Generic matrices",
    "title": "Inverse",
    "category": "section",
    "text": "inv{T <: RingElem}(::MatElem{T})\ninv{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(QQ, \"x\")\nK, a = NumberField(x^3 + 3x + 1, \"a\")\nS = MatrixSpace(K, 3, 3)\n\nA = S([K(0) 2a + 3 a^2 + 1; a^2 - 2 a - 1 2a; a^2 + 3a + 1 2a K(1)])\n\nX = inv(A)\n\nR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 3, 3)\n\nA = S([R(0) 2x + 3 x^2 + 1; x^2 - 2 x - 1 2x; x^2 + 3x + 1 2x R(1)])\n    \nX, d = inv(A)"
},

{
    "location": "matrix/#LinearAlgebra.nullspace-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.nullspace",
    "category": "method",
    "text": "nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: RingElement}\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In AbstractAlgebra we use the name ``kernel\'\' for a function to compute an integral kernel.\n\n\n\n"
},

{
    "location": "matrix/#LinearAlgebra.nullspace-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.nullspace",
    "category": "method",
    "text": "nullspace(M::AbstractAlgebra.MatElem{T}) where {T <: RingElement}\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In AbstractAlgebra we use the name ``kernel\'\' for a function to compute an integral kernel.\n\n\n\nnullspace(M::AbstractAlgebra.MatElem{T}) where {T <: FieldElement}\n\nReturns a tuple (nu N) consisting of the nullity nu of M and a basis N (consisting of column vectors) for the right nullspace of M, i.e. such that MN is the zero matrix. If M is an mtimes n matrix N will be an ntimes nu matrix. Note that the nullspace is taken to be the vector space kernel over the fraction field of the base ring if the latter is not a field. In Nemo we use the name ``kernel\'\' for a function to compute an integral kernel.\n\n\n\n"
},

{
    "location": "matrix/#Nullspace-1",
    "page": "Generic matrices",
    "title": "Nullspace",
    "category": "section",
    "text": "nullspace{T <: RingElem}(::MatElem{T})\nnullspace{T <: FieldElem}(::MatElem{T})ExamplesR, x = PolynomialRing(ZZ, \"x\")\nS = MatrixSpace(R, 4, 4)\n   \nM = S([-6*x^2+6*x+12 -12*x^2-21*x-15 -15*x^2+21*x+33 -21*x^2-9*x-9;\n       -8*x^2+8*x+16 -16*x^2+38*x-20 90*x^2-82*x-44 60*x^2+54*x-34;\n       -4*x^2+4*x+8 -8*x^2+13*x-10 35*x^2-31*x-14 22*x^2+21*x-15;\n       -10*x^2+10*x+20 -20*x^2+70*x-25 150*x^2-140*x-85 105*x^2+90*x-50])\n   \nn, N = nullspace(M)"
},

{
    "location": "matrix/#LinearAlgebra.hessenberg-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "LinearAlgebra.hessenberg",
    "category": "method",
    "text": "hessenberg(A::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturns the Hessenberg form of M, i.e. an upper Hessenberg matrix which is similar to M. The upper Hessenberg form has nonzero entries above and on the diagonal and in the diagonal line immediately below the diagonal.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.ishessenberg-Union{Tuple{MatElem{T}}, Tuple{T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.ishessenberg",
    "category": "method",
    "text": "ishessenberg(A::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturns true if M is in Hessenberg form, otherwise returns false.\n\n\n\n"
},

{
    "location": "matrix/#Hessenberg-form-1",
    "page": "Generic matrices",
    "title": "Hessenberg form",
    "category": "section",
    "text": "hessenberg{T <: RingElem}(::MatElem{T})ishessenberg{T <: RingElem}(::MatElem{T})ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\n   \nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nA = hessenberg(M)\nishessenberg(A) == true"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.charpoly-Union{Tuple{T}, Tuple{Ring,MatElem{T}}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.charpoly",
    "category": "method",
    "text": "charpoly(V::Ring, Y::Generic.MatrixElem{T}) where {T <: RingElement}\n\nReturns the characteristic polynomial p of the matrix M. The polynomial ring R of the resulting polynomial must be supplied and the matrix is assumed to be square.\n\n\n\n"
},

{
    "location": "matrix/#Characteristic-polynomial-1",
    "page": "Generic matrices",
    "title": "Characteristic polynomial",
    "category": "section",
    "text": "charpoly{T <: RingElem}(::Ring, ::MatElem{T})ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\nT, x = PolynomialRing(R, \"x\")\n\nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nA = charpoly(T, M)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.minpoly-Union{Tuple{T}, Tuple{Ring,MatElem{T},Bool}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.minpoly",
    "category": "method",
    "text": "minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring S of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.minpoly-Union{Tuple{T}, Tuple{Ring,MatElem{T},Bool}} where T<:FieldElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.minpoly",
    "category": "method",
    "text": "minpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: FieldElement}\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring S of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\nminpoly(S::Ring, M::MatElem{T}, charpoly_only::Bool = false) where {T <: RingElement}\n\nReturns the minimal polynomial p of the matrix M. The polynomial ring S of the resulting polynomial must be supplied and the matrix must be square.\n\n\n\n"
},

{
    "location": "matrix/#Minimal-polynomial-1",
    "page": "Generic matrices",
    "title": "Minimal polynomial",
    "category": "section",
    "text": "minpoly{T <: RingElem}(::Ring, ::MatElem{T}, ::Bool)\nminpoly{T <: FieldElem}(::Ring, ::MatElem{T}, ::Bool)ExamplesR = GF(13)\nT, y = PolynomialRing(R, \"y\")\n   \nM = R[7 6 1;\n      7 7 5;\n      8 12 5]\n\nA = minpoly(T, M)"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.similarity!-Union{Tuple{T}, Tuple{MatElem{T},Int64,T}} where T<:RingElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.similarity!",
    "category": "method",
    "text": "similarity!(A::Generic.MatrixElem{T}, r::Int, d::T) where {T <: RingElement}\n\nApplies a similarity transform to the ntimes n matrix M in-place. Let P be the ntimes n identity matrix that has had all zero entries of row r replaced with d, then the transform applied is equivalent to M = P^-1MP. We require M to be a square matrix. A similarity transform preserves the minimal and characteristic polynomials of a matrix.\n\n\n\n"
},

{
    "location": "matrix/#Transforms-1",
    "page": "Generic matrices",
    "title": "Transforms",
    "category": "section",
    "text": "similarity!{T <: RingElem}(::MatElem{T}, ::Int, ::T)ExamplesR = ResidueRing(ZZ, 7)\nS = MatrixSpace(R, 4, 4)\n   \nM = S([R(1) R(2) R(4) R(3); R(2) R(5) R(1) R(0);\n       R(6) R(1) R(3) R(2); R(1) R(1) R(3) R(5)])\n   \nsimilarity!(M, 1, R(3))"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.weak_popov-Union{Tuple{Mat{T}}, Tuple{T}} where T<:PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.weak_popov",
    "category": "method",
    "text": "weak_popov(A::Mat{T}) where {T <: PolyElem}\n\nReturn the weak Popov form of A.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.weak_popov_with_trafo-Union{Tuple{Mat{T}}, Tuple{T}} where T<:PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.weak_popov_with_trafo",
    "category": "method",
    "text": "weak_popov_with_trafo(A::Mat{T}) where {T <: PolyElem}\n\nCompute a tuple (P U) where P is the weak Popov form of A and U is a transformation matrix so that P = UA.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.popov-Union{Tuple{Mat{T}}, Tuple{T}} where T<:PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.popov",
    "category": "method",
    "text": "popov(A::Mat{T}) where {T <: PolyElem}\n\nReturn the Popov form of A.\n\n\n\n"
},

{
    "location": "matrix/#AbstractAlgebra.Generic.popov_with_trafo-Union{Tuple{Mat{T}}, Tuple{T}} where T<:PolyElem",
    "page": "Generic matrices",
    "title": "AbstractAlgebra.Generic.popov_with_trafo",
    "category": "method",
    "text": "popov_with_trafo(A::Mat{T}) where {T <: PolyElem}\n\nCompute a tuple (P U) where P is the Popov form of A and U is a transformation matrix so that P = UA.\n\n\n\n"
},

{
    "location": "matrix/#(Weak)-Popov-form-1",
    "page": "Generic matrices",
    "title": "(Weak) Popov form",
    "category": "section",
    "text": "AbstractAlgebra.jl provides algorithms for computing the (weak) Popov of a matrix with entries in a univariate polynomial ring over a field.weak_popov{T <: PolyElem}(::Generic.Mat{T})\nweak_popov_with_trafo{T <: PolyElem}(::Generic.Mat{T})\npopov{T <: PolyElem}(::Generic.Mat{T})\npopov_with_trafo{T <: PolyElem}(::Generic.Mat{T})"
},

{
    "location": "matrix_algebras/#",
    "page": "Generic matrix algebras",
    "title": "Generic matrix algebras",
    "category": "page",
    "text": "CurrentModule = AbstractAlgebra"
},

{
    "location": "matrix_algebras/#Generic-matrix-algebras-1",
    "page": "Generic matrix algebras",
    "title": "Generic matrix algebras",
    "category": "section",
    "text": "AbstractAlgebra.jl allows the creation of an algebra (ring) of mtimes m matrices over a computable, commutative ring.Functions specific to generic matrix algebras of mtimes m matrices are implemented in  src/generic/MatrixAlgebra.jl. The remaining functionality is in the file src/generic/Matrix.jl.As well as implementing the entire Matrix interface, including the optional functionality, there are many additional generic algorithms implemented for matrix algebras.Almost all of the functionality specified for generic matrices is available for matrix algebras. The exceptions are functions such as solve and nullspace which may return non-square matrices, or which don\'t accept square matrices.All of the generic functionality is part of the Generic submodule of AbstractAlgebra.jl. This is exported by default, so it is not necessary to qualify names of functions."
},

{
    "location": "matrix_algebras/#Types-and-parent-objects-1",
    "page": "Generic matrix algebras",
    "title": "Types and parent objects",
    "category": "section",
    "text": "Generic matrices in AbstractAlgebra.jl have type Generic.MatAlgElem{T} for matrices in  a matrix algebra, where T is the type of elements of the matrix. Internally, generic matrices are implemented using an object wrapping a Julia two dimensional array, though they are not themselves Julia arrays. See the file src/generic/GenericTypes.jl for details.Parents of generic matrices in a matrix algebra have type Generic.MatAlgebra{T}.Note that matrix algebras are noncommutative rings. Thus their types belong to NCRing and NCRingElem. They cannot be used in constructions which require a commutative ring (Ring and RingElem respectively).The generic matrix algebra matrix types belong to the abstract type AbstractAlgebra.MatAlgElem{T} and the parent types belong to  AbstractAlgebra.MatAlgebra{T} Note that both of these require disambiguation from the concrete types in Generic of the same name.The degree and base ring R of a generic matrix are stored in its parent object, however to allow creation of matrices without first creating the matrix space parent, generic matrices in Julia do not contain a reference to their parent. They contain the row and column numbers (or degree, in the case of matrix algebras) and the base ring on a per matrix basis. The parent object can then be reconstructed from this data on demand."
},

{
    "location": "matrix_algebras/#Matrix-space-constructors-1",
    "page": "Generic matrix algebras",
    "title": "Matrix space constructors",
    "category": "section",
    "text": "A matrix algebra in AbstractAlgebra.jl represents a collection of all matrices with given degree and base ring.In order to construct matrices in AbstractAlgebra.jl, one must construct the matrix algebra itself. This is accomplished with the following constructor.MatrixAlgebra(R::Ring, degree::Int; cache::Bool=true)Construct the algebra of matrices with the given degree over the given base ring. By default such matrix spaces are cached based on the base ring and degree. If the optional named parameter cached is set to false, no caching occurs.Here are some examples of creating matrix algebras and making use of the resulting parent objects to coerce various elements into the matrix algebra.ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixAlgebra(R, 3)\n\nA = S()\nB = S(12)\nC = S(R(11))"
},

{
    "location": "matrix_algebras/#Matrix-algebra-functionality-provided-by-AbstractAlgebra.jl-1",
    "page": "Generic matrix algebras",
    "title": "Matrix algebra functionality provided by AbstractAlgebra.jl",
    "category": "section",
    "text": "Most of the generic matrix functionality described in the generic matrix section of the documentation is available for both matrix spaces and matrix algebras. Exceptions include functions that do not return or accept square matrices or which cannot specify a parent. Such functions include solve and nullspace which can\'t be provided for matrix algebras.In addition to the functionality described for matrix spaces, matrix algebras support all noncommutative ring operations, and matrix algebras can be used as a base ring for other generic constructs that accept a noncommutative base ring (NCRing).In this section we describe functionality provided for matrix algebras only."
},

{
    "location": "matrix_algebras/#AbstractAlgebra.Generic.degree-Tuple{AbstractAlgebra.Generic.MatAlgElem}",
    "page": "Generic matrix algebras",
    "title": "AbstractAlgebra.Generic.degree",
    "category": "method",
    "text": "degree(a::Generic.MatAlgElem)\n\nReturn the degree n of the ntimes n matrix a..\n\n\n\n"
},

{
    "location": "matrix_algebras/#Basic-matrix-functionality-1",
    "page": "Generic matrix algebras",
    "title": "Basic matrix functionality",
    "category": "section",
    "text": "As well as the Ring and Matrix interfaces, the following functions are provided to manipulate matrices.degree(::Generic.MatAlgElem)ExamplesR, t = PolynomialRing(QQ, \"t\")\nS = MatrixAlgebra(R, 3)\n\nA = S([t + 1 t R(1); t^2 t t; R(-2) t + 2 t^2 + t + 1])\n\nn = degree(A)"
},

{
    "location": "map/#",
    "page": "Map interface",
    "title": "Map interface",
    "category": "page",
    "text": ""
},

{
    "location": "map/#Map-interface-1",
    "page": "Map interface",
    "title": "Map interface",
    "category": "section",
    "text": "Maps in Nemo model maps on sets f  D to C for some domain D and codomain C, which have no real limitations except that elements of the codomain and domain be represented by element objects in the system..Maps f  D to C in AbstractAlgebra are modeled by Julia objects that are able to be called on a single element c in C of the domain to yield an element f(c) of the codomain. We say that the map is being applied.Maps can be constructed from Julia functions, or they can be represented by some other kind of data, e.g. a matrix, or built up from other maps.In the following, we will always use the word \"function\" to mean a Julia function, and reserve the word \"map\" for a map on sets, whether mathematically, or as an object in the system..Maps in AbstractAlgebra have a domain and codomain, can be applied, composed and composed with the identity map (assuming its domain is compatible). Various special kinds of map provide more functionality."
},

{
    "location": "map/#Parent-objects-1",
    "page": "Map interface",
    "title": "Parent objects",
    "category": "section",
    "text": "Maps in AbstractAlgebra currently don\'t have parents. This will change later when AbstractAlgebra has a category system, so that the parent of a map can be some sort of Hom set."
},

{
    "location": "map/#Map-classes-1",
    "page": "Map interface",
    "title": "Map classes",
    "category": "section",
    "text": "All maps in AbstractAlgebra belong to a class of maps. The classes are modeled as abstract types that lie in a hierarchy, inheriting from SetMap at the top of the hierarchy. Other classes that inherit from SetMap are FunctionalMap for maps that are constructed from a Julia function (or closure), and IdentityMap for the class of the identity maps within the system.One might naturally assume that map types belong directly to these classes in the way that types of other objects in the system belong to abstract types in the AbstractAlgebra type hierarchy. However, in order to provide an extensible system, this is not the case.Instead, a map type MyMap will belong to an abstract type of the form Map{D, C, T, MyMap}, where D is the type of the object representing the domain of the map type (this can also be an abstract type, such as Group), C is the type of the object representing the codomain of the map type and T is the map class that MyMap belongs to, e.g. SetMap or FunctionalMap.Because a four parameter type system becomes quite cumbersome to use, we provide a number of functions for referring to collections of map types.If writing a function that accepts any map type, one makes the type of its argument belong to Map. For example f(M::Map) = 1.If writing a function that accepts any map from a domain of type D to a codomain of type C, one makes writes for example f(M::Map{D, C}) = 2. Note that D and C can be abstract types, such as Group, but otherwise must be the types of the parent objects representing the domain and codomain.A function that accepts any map belonging to a given map class might be written as f(M::Map(FunctionalMap)) = 3 or f(M::Map(FunctionalMap){D, C}) = 4 for example, where D and C are the types of the parent objects for the domain and codomain.Finally, if a function should only work for a map of a given map type MyMap, say, one writes this f(M::Map(MyMap)) or f(M::Map(MyMap){D, C}, where as usual D and C are the types of the domain and codomain parent objects."
},

{
    "location": "map/#Implementing-new-map-types-1",
    "page": "Map interface",
    "title": "Implementing new map types",
    "category": "section",
    "text": "There are two common kinds of map type that developers will need to write. The first has a fixed domain and codomain, and the second is a type parameterised by the types of the domain and codomain. We give two simple examples here of how this might look.In the case of fixed domain and codomain, e.g. Integers{BigInt}, we would write it as follows:mutable struct MyMap <: Map{Integers{BigInt}, Integers{BigInt}, SetMap, MyMap}\n   # some data fields\nendIn the case of parameterisation by the type of the domain and codomain:mutable struct MyMap{D, C} <: Map{D, C, SetMap, MyMap}\n   # some data fields\nendAs mentioned above, to write a function that only accepts maps of type MyMap, one writes the functions as follows:function my_fun(M::Map(MyMap))The Map function then computes the correct type to use, which is actually not MyMap if all features of the generic Map infrastructure are required. It is bad practice to write functions for MyMap directly instead of Map(MyMap), since other users will be unable to use generic constructions over the map type MyMap."
},

{
    "location": "map/#Getters-and-setters-1",
    "page": "Map interface",
    "title": "Getters and setters",
    "category": "section",
    "text": "When writing new map types, it is very important to define getters and setters of the fields of the new map type, rather than to access them directly.Let us suppose that the MyMap type has a field called foo. Rather than access this field by writing M.foo, one must access it using foo(M) (at least until Julia 1.1).If such a getter only needs to access the field foo of M, there is a standard way of defining such a getter and setter when defining a new map type.foo(M::Map(MyMap)) = get_field(M, :foo)To set a field of a map, one needs a setter, which can be implemented as follows:set_foo!(M::Map(MyMap), a) = set_field(M, :foo, a)In general, setters should be used rarely for map types.Note that by providing getter and setter functions, map types need not even contain fields with the given name. For example, for a MyMap map type for maps between integers, one does not wish to explicitly store the domain and codomain in MyMap. Instead, we can define the getter functions domain and codomain to return JuliaZZ for any MyMap object. domain(M::Map(MyMap)) = JuliaZZ\ncodomain(M::Map(MyMap)) = JuliaZZ"
},

{
    "location": "map/#Required-functionality-for-maps-1",
    "page": "Map interface",
    "title": "Required functionality for maps",
    "category": "section",
    "text": "All map types must implement a standard interface, which we specify here.We will define this interface for a custom map type MyMap belonging to Map(SetMap), SetMap being the map class that all maps types belong to.Note that map types do not need to contain any specific fields, but must provide accessor functions (getters and setters) in the manner described above.The required accessors for map types of class SetMap are as follows.domain(M::Map(MyMap))\ncodomain(M::Map(MyMap))Return the domain and codomain parent objects respectively, for the map M. it is only necessary to define these functions if the map type MyMap does not contain fields domain and codomain containing these parent objects.It is also necessary to be able to apply a map. This amounts to overloading the call method for objects belonging to Map(MyMap).(M::Map(MyMap)(a))Apply the map M to the element a of the domain of M. Note that it is usual to add a type assertion to the return value of this function, asserting that the return value has type elem_type(C) where C is the type of the codomain parent object."
},

{
    "location": "map/#Optional-functionality-for-maps-1",
    "page": "Map interface",
    "title": "Optional functionality for maps",
    "category": "section",
    "text": "The Generic module in AbstractAlgebra automatically provides certain functionality for map types, assuming that they satisfy the full interface described above.However, certain map types or map classes might like to provide their own implementation of this functionality, overriding the generic functionality.We describe this optional functionality here."
},

{
    "location": "map/#Show-method-1",
    "page": "Map interface",
    "title": "Show method",
    "category": "section",
    "text": "Custom map types may like to provide a custom show method if the default of displaying the domain and codomain of the map is not sufficient.show(io::IO, M::Map(MyMap))"
},

{
    "location": "map/#Identity-maps-1",
    "page": "Map interface",
    "title": "Identity maps",
    "category": "section",
    "text": "There is a concrete map type Generic.IdentityMap{D} for the identity map on a given domain. Here D is the type of the object representing that domain.Generic.IdentityMap belongs to the supertype Map{D, C, AbstractAlgebra.IdentityMap, IdentityMap}.Note that the map class is also called IdentityMap. It is an abstract type, whereas Generic.IdentityMap is a concrete type in the Generic module.An identity map has the property that when composed with any map whose domain or codomain is compatible, that map will be returned as the composition. Identity maps can therefore serve as a starting point when building up a composition of maps, starting an identity map.We do not cached identity maps in the system, so that if more than one is created on the same domain, there will be more than one such map in the system. This underscores the fact that there is in general no way for the system to know if two maps compose to give an identity map, and therefore the only two maps that can be composed to give an identity map are identity maps on the same domain.To construct an identity map for a given domain, specified by a parent object R, say, we have the following function.identity_map(R::D) where D <: SetReturn an identity map on the domain R.Of course there is nothing stopping a map type or class from implementing its own identity map type, and defining composition of maps of the same kind with such an identity map. In such a case, the class of such an identity map type must belong to AbstractAlgebra.IdentityMap so that composition with other map types still works."
},

{
    "location": "map/#Composition-of-maps-1",
    "page": "Map interface",
    "title": "Composition of maps",
    "category": "section",
    "text": "Any two compatible maps in AbstractAlgebra can be composed and any composition can be applied.In order to facilitate this, the Generic module provides a type CompositionMap{D, C}, which contains two maps map1 and map2, corresponding to the two maps to be applied in a composition, in the order they should be applied.To construct a composition map from two existing maps, we have the following function:compose(f::Map{D, U}, g::Map{U, C}) where {D, U, C}Compose the two maps f and g, i.e. return the map h such that h(x) = g(f(x)).As a shortcut for this function we have the following operator:*(f::Map{D, U}, g::Map{U, C}) where {D, U, C} = compose(f, g)Note the order of composition. If we have maps f  X to Y, g  Y to Z the correct order of the maps in this operator is f*g, so that (f*g)(x) = g(f(x)).This is chosen so that for left R-module morphisms represented by a matrix, the order of matrix multiplication will match the order of composition of the corresponding morphisms.Of course, a custom map type or class of maps can implement its own composition type and compose function.This is the case with the FunctionalMap class for example, which caches the Julia function/closure corresponding to the composition of two functional maps. As this cached function needs to be stored inside the composition, a special type is necessary for the composition of two functional maps.By default, compose will check that the two maps are composable, i.e. the codomain of the first map matches the domain of the second map. This is implemented by the following function:check_composable(f::Map{D, U}, g::Map{U, C})Raise an exception if the codomain of f doesn\'t match the domain of g.Note that composite maps should keep track of the two maps they were constructed from. To access these maps, the following functions are provided:map1(f::CompositeMap)\nmap2(f::CompositeMap)Any custom composite map type must also provide these functions for that map type, even if there exist fields with those names. This is because there is no common map class for all composite map types. Therefore the Generic system cannot provide fallbacks for all such composite map types."
},

{
    "location": "functional_map/#",
    "page": "Functional maps",
    "title": "Functional maps",
    "category": "page",
    "text": ""
},

{
    "location": "functional_map/#Functional-maps-1",
    "page": "Functional maps",
    "title": "Functional maps",
    "category": "section",
    "text": "A functional map in AbstractAlgebra is a map which can be applied by evaluating a Julia function or closure. It is represented by a map object that contains such a function/closure, usually in a field called image_fn.All functional maps belong to the map class AbstractAlgebra.FunctionalMap.A generic concrete type Generic.FunctionalMap is provided by the Generic module to implement a generic functional map type. This allows for functional maps that contain no extra data, other than a Julia function/closure.Custom map types can also be defined which have map class AbstractAlgebra.FunctionalMap."
},

{
    "location": "functional_map/#Functional-map-interface-1",
    "page": "Functional maps",
    "title": "Functional map interface",
    "category": "section",
    "text": "All functional map types must define their supertypes as in the following example:mutable struct MyFunctionalMap{D, C} <: Map{D, C, FunctionalMap, MyFunctionalMap}\n   # some fields\n   image_fn::Function\nendOf course MyFunctionalMap need not be parameterised if the types D and C of the domain and codomain objects are known."
},

{
    "location": "functional_map/#Required-functions-for-functional-maps-1",
    "page": "Functional maps",
    "title": "Required functions for functional maps",
    "category": "section",
    "text": "The following functions must be defined for all functional map types or classes:image_fn(M::Map(MyFunctionalMap))Return the Julia function or closure that corresponds to application of the map M. This function only needs to be provided if this function is not stored in an image_fn field of the MyFunctionalMap type."
},

{
    "location": "functional_map/#Generic-functional-maps-1",
    "page": "Functional maps",
    "title": "Generic functional maps",
    "category": "section",
    "text": "The Generic module provides a concrete type FunctionalMap which merely keeps track of a Julia function/closure implementing the map.Such maps can be constructed using the following function:map_from_func(R, S, f::Function)Construct the generic functional map with domain and codomain given by the parent objects R and S corresponding to the Julia function f.Examplesf = map_from_func(ZZ, ZZ, x -> x + 1)\n\nf(ZZ(2))"
},

{
    "location": "map_cache/#",
    "page": "Cached maps",
    "title": "Cached maps",
    "category": "page",
    "text": ""
},

{
    "location": "map_cache/#Cached-maps-1",
    "page": "Cached maps",
    "title": "Cached maps",
    "category": "section",
    "text": "All basic map (i.e. those not built up from other maps) in AbstractAlgebra can be cached.A cache is a dictionary that can be switched on and off at run time that keeps a cache of previous evaluations of the map. This can be useful if the map is extremely difficult to evaluate, e.g. a discrete logarithm map. Rather than evaluate the map afresh each time, the map first looks up the dictionary of previous known values of the map.To facilitate caching of maps, the Generic module provides a type Generic.MapCache, which can be used to wrap any existing map object with a dictionary.Importantly, the supertype of the resulting MapCache object is identical to that of the map being cached. This means that any functions that would accept the original map will also accept the cached version.Note that caching of maps only works for maps that correctly abstract access to their fields using accessor functions, as described in the map interface."
},

{
    "location": "map_cache/#Cached-map-constructors-1",
    "page": "Cached maps",
    "title": "Cached map constructors",
    "category": "section",
    "text": "To construct a cached map from an existing map object, we have the following function:cached(M::Map; enabled=true, limit=100)Return a cached map with the same supertype as M, caching up to limit values of the map M in a dictionary, assuming that the cache is enabled.Caches can be disabled by setting the value of the parameter enabled to false. This allows for the user to quickly go through code and completely disable caches of maps that were previously enabled, for testing purposed, etc.Caches can also be turned on and off at run time (see below).Examplesf = map_from_func(ZZ, ZZ, x -> x + 1)\ng = cached(f)\n\nf(ZZ(1)) == g(ZZ(1))"
},

{
    "location": "map_cache/#Functionality-for-cached-maps-1",
    "page": "Cached maps",
    "title": "Functionality for cached maps",
    "category": "section",
    "text": "The following functions are provided for cached maps.enable_cache!(M::MapCache)\ndisable_cache!(M::MapCache)Temporarily enable or disable the cache for the given map. The values stored in the cache are not lost when it is disabled.set_limit!(M::MapCache, limit::Int)Set the limit on the number of values that can be cached in the dictionary, to the given value. Setting the value to 0 will effectively disable further caching for this map.Examplesf = cached(map_from_func(ZZ, ZZ, x -> x + 1))\n\na = f(ZZ(1))\ndisable_cache!(f)\nb = f(ZZ(1))\nenable_cache!(f)\nc = f(ZZ(1))\n\nset_limit!(f, 200)\nd = f(ZZ(1))"
},

{
    "location": "map_with_inverse/#",
    "page": "Map with inverse",
    "title": "Map with inverse",
    "category": "page",
    "text": ""
},

{
    "location": "map_with_inverse/#Map-with-inverse-1",
    "page": "Map with inverse",
    "title": "Map with inverse",
    "category": "section",
    "text": "It is not possible to provide generic functionality to invert a map. However, sometimes one knows an inverse map explicitly and would like to keep track of this.Recall that as map composition is not commutative, there is a notion of a left inverse and a right inverse for maps.To keep track of such inverse maps, AbstractAlgebra provides data types Generic.MapWithRetraction and GenericMapWithSection.Given a map f  X to Y, a retraction of f is a map g  Y to X such that g(f(x)) = x for all x in X.Given a map f  X to Y, a section of f is a map g  Y to X such that f(g(x)) = x for all y in Y.In AbstractAlgebra, a map with retraction/section is an object containing a pair of maps, the second of which is a retraction/section of the first.Maps with retraction/section can be composed, and we also define the inverse of such a pair to be the map with the pair swapped. Thus the inverse of a map with retraction is a map with section. "
},

{
    "location": "map_with_inverse/#Map-with-inverse-constructors-1",
    "page": "Map with inverse",
    "title": "Map with inverse constructors",
    "category": "section",
    "text": "To construct a map with retraction/section from a pair of maps, we have the following functions:map_with_retraction(m::Map{D, C}, r::Map{C, D}) where {D, C}\nmap_with_section(m::Map{D, C}, s::Map{C, D}) where {D, C}Construct the map with retraction/section given a known retraction/section r or s respectively, of m.For convenience we allow construction of maps with retraction/section from a pair of Julia functions/closures.map_with_retraction_from_func(R, S, f::Function, r::Function)\nmap_with_section_from_func(R, S, f::Function, s::Function)Construct the map with retraction/section such that the map is given by the function f and the retraction/section is given by the function r or s respectively. Here R is the parent object representing the domain and S is the parent object representing the codomain of f.Examplesf = map_with_retraction_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)\n\na = f(ZZ(1))"
},

{
    "location": "map_with_inverse/#Functionality-for-maps-with-inverses-1",
    "page": "Map with inverse",
    "title": "Functionality for maps with inverses",
    "category": "section",
    "text": "The following functionality is provided for maps with inverses.inv(M::MapWithRetraction)\ninv(M::MapWithSection)Return the map with the two maps contained in M swapped. In the first case, a MapWithSection is returned. In the second case a MapWithRetraction is returned.To access the two maps stored in a map with retraction/section, we have the following:image_map(M::MapWithRetraction)\nimage_map(M::MapWithSection)\nretraction_map(M::MapWithRetraction)\nsection_map(M::MapWithSection)The first two of these functions return the first map in a map with retraction/section, the second two functions return the corresponding second maps.Examplesf = map_with_retraction_from_func(ZZ, ZZ, x -> x + 1, x -> x - 1)\ng = inv(f)\nh = f*g\n\na = h(ZZ(1))"
},

{
    "location": "types/#",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Appendix A: Types in AbstractAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "types/#Appendix-A:-Types-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Appendix A: Types in AbstractAlgebra.jl",
    "category": "section",
    "text": "On this page we discuss the abstract type hierarchy in AbstractAlgebra.jl and objects known as parents which contain additional information about groups, rings, fields and modules, etc., that can\'t be stored in types alone.These details are technical and can be skipped or skimmed by new users of  Julia/AbstractAlgebra.jl. Types are almost never dealt with directly when scripting  AbstractAlgebra.jl to do mathematical computations. In contrast, AbstractAlgebra.jl developers will want to know how we model mathematical objects and their rings, fields, groups, etc."
},

{
    "location": "types/#The-abstract-type-hierarchy-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "The abstract type hierarchy in AbstractAlgebra.jl",
    "category": "section",
    "text": "Abstract types in Julia can also belong to one another in a hierarchy. We make use of such a hierarchy to organise the kinds of mathematical objects in AbstractAlgebra.jl.For example, the AbstractAlgebra.Field abstract type belongs to the  AbstractAlgebra.Ring abstract type. In practice, this means that any generic function in AbstractAlgebra.jl which is designed to work with ring objects will also work with field objects.In AbstractAlgebra.jl we also distinguish between the elements of a field, say, and the field itself.For example, we have an object of type Generic.PolyRing to model a generic polynomial ring, and elements of that polynomial ring would have type Generic.Poly. For this purpose, we also have a hierarchy of abstract types, such as FieldElem, that the types of element objects can belong to.(Image: alt text)"
},

{
    "location": "types/#Why-types-aren\'t-enough-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Why types aren\'t enough",
    "category": "section",
    "text": "Naively, one might have expected that rings in AbstractAlgebra.jl could be modeled as types and their elements as objects with the given type. But there are various reasons why this is not a good model.Consider the ring R = mathbbZnmathbbZ for a multiprecision integer n. If we were to model the ring R as a type, then the type would somehow need to contain the modulus n. This is not possible in Julia, and in fact it is not desirable, since the compiler would then recompile all the associated functions every time a different modulus n was used.We could attach the modulus n to the objects representing elements of the ring, rather than their type.But now we cannot create new elements of the ring mathbbZnmathbbZ given only their type, since the type no longer contains the modulus n.Instead, the way we get around this in AbstractAlgebra.jl is to have special (singleton) objects that act like types, but are really just ordinary Julia objects. These objects, called parent objects can contain extra information, such as the modulus n. In order to create new elements of mathbbZnmathbbZ as above, we overload the call operator for the parent object.In the following AbstractAlgebra.jl example, we create the parent object R corresponding to the ring mathbbZ7mathbbZ. We then create a new element a of this ring by calling the parent object R.R = ResidueRing(ZZ, 7)\na = R(3)Here, R is the parent object, containing the modulus 7. So this example creates  the element a = 3 pmod7."
},

{
    "location": "types/#More-complex-example-of-parent-objects-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "More complex example of parent objects",
    "category": "section",
    "text": "Here is some Julia/AbstractAlgebra.jl code which constructs a polynomial ring over the integers, a polynomial in that ring and then does some introspection to illustrate the various relations between the objects and types.using AbstractAlgebra\n\nR, x = ZZ[\"x\"]\n\nf = x^2 + 3x + 1\n\ntypeof(R) <: PolyRing\n\ntypeof(f) <: PolyElem\n\nparent(f) == R"
},

{
    "location": "types/#Concrete-types-in-AbstractAlgebra.jl-1",
    "page": "Appendix A: Types in AbstractAlgebra.jl",
    "title": "Concrete types in AbstractAlgebra.jl",
    "category": "section",
    "text": "Here we give a list of the concrete types in AbstractAlgebra.jl.In parentheses we put the types of the corresponding parent objects.perm{<:Integer} (PermGroup{<:Integer})\ngfelem{<:Integer} (GFField{<:Integer})We also think of various Julia types as though they were AbstractAlgebra.jl types:BigInt (Integers{BigInt})\nRational{BigInt} (Rationals{BigInt})Then there are various types for generic constructions over a base ring. They are all parameterised by a type T which is the type of the elements of the base ring they are defined over. Generic.Poly{T} (Generic.PolyRing{T})\nGeneric.MPoly{T} (Generic.MPolyRing{T})\nGeneric.RelSeries{T} (Generic.RelSeriesRing{T})\nGeneric.AbsSeries{T} (Generic.AbsSeriesRing{T})\nGeneric.LaurentSeriesRingElem{T} (Generic.LaurentSeriesRing{T})\nGeneric.LaurentSeriesFieldElem{T} (Generic.LaurentSeriesField{T})\nGeneric.Res{T} (Generic.ResRing{T})\nGeneric.Frac{T} (Generic.FracField{T})\nGeneric.Mat{T} (Generic.MatSpace{T})"
},

]}
