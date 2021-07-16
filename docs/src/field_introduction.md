# Introduction.

A number of basic fields are provided, such as the rationals, finite
fields, the real field, etc.

Various generic field constructions can then be made recursively on
top of these basic fields. For example, fraction fields, residue fields,
function fields, etc.

From the point of view of the system, all fields are rings and whether
an object is a ring/field or an element thereof can be determined at
the type level. There are abstract types for all field and for all
field element types.

The field hierarchy can be extended by implementing new fields to follow
one or more field interfaces, including the interface that all fields
must follow. Once an interface is satisfied, all the corresponding
generic functionality will work over the new field.

Implementations of new fields can either be generic or can be specialised
implementations provided by, for example, a C library.



