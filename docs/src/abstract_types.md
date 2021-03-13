# Abstract types

AbstractAlgebra.jl implements a couple of abstract types which can be extended.

## Parents

```@setup parents
using Kroki

parents_bad = plantuml"""
@startwbs
skinparam monochrome true
skinparam arrowThickness 3
skinparam WbsType tree
<style>
wbsDiagram {
  FontName Monospaced
  FontSize 16
  RoundCorner 10
  Padding 5
  Margin 15
}
</style>

+ Set
++ Ideal{T}
++ Group
--- AbstractPermutationGroup
++ AdditiveGroup
--- Module{T}
++++ MatSpace{T}
---- FPModule{T}
++ NCRing
+++ NCPolyRing{T}
--- MatAlgebra{T}
--- Ring
++++ PolyRing{T}
++++ LaurentPolynomialRing{T}
++++ SeriesRing{T}
++++ ResRing{T}
++++ Field
+++++ ResField{T}
+++++ FracField{T}
+++++ NumField{T}
++++++ SimpleNumField{T}
+++++ FinField
@endwbs
"""

parents = plantuml"""
@startuml
skinparam monochrome true
skinparam defaultFontName Monospaced
skinparam defaultFontSize 16
skinparam objectArrowColor DarkGray
skinparam RoundCorner 15

' Ideals
"Set" -down-> "Ideal{T}"

' Groups
"Set" -down-> "Group"
"Group" -down-> "AbstractPermutationGroup"

' Additive groups
"Set" -down-> "AdditiveGroup"
"AdditiveGroup" -down-> "Module{T}"
"Module{T}" -down-> "FPModule{T}"
"Module{T}" -down-> "MatSpace{T}"

' Non-commutative rings
"Set" -down-> "NCRing"
"NCRing" -down--> "MatAlgebra{T}"
"NCRing" -down--> "NCPolyRing{T}"
"NCRing" -down-> "Ring"

' Rings
"Ring" -down-> "PolyRing{T}"
"PolyRing{T}" -down-> "LaurentPolynomialRing{T}"
"PolyRing{T}" -down-> "SeriesRing{T}"
"Ring" -down-> "ResRing{T}"
"Ring" -down---> "Field"

' Fields
"Field" -down-> "ResField{T}"
"Field" -down-> "FracField{T}"
"Field" -down-> "NumField{T}"
"NumField{T}" -down-> "SimpleNumField{T}"
"Field" -down-> "FinField"

hide members
hide circle

@enduml
"""
```
```@example parents
parents #hide
```

## Elements

```@setup elements
using Kroki

elements = plantuml"""
@startuml
skinparam monochrome true
skinparam defaultFontName Monospaced
skinparam defaultFontSize 16
skinparam objectArrowColor DarkGray
skinparam RoundCorner 15


' Ideals
"SetElem" -down-> "IdealElem{T}"

' Maps
' NOTE: Identity is 
"SetElem" -down----> "Map{D, C, S, T}"
"Map{D, C, S, T}" .right. "SetMap"
"SetMap" -down-> "IdentityMap"
"SetMap" -down-> "FunctionalMap"
"FunctionalMap" -down-> "FPModuleHomomorphism"

' Groups
' NOTE: Shouldn't AbstractPerm have the suffix -Elem?
"SetElem" -down-> "GroupElem"
"GroupElem" -down-> "AbstractPerm"

' Additive groups
"SetElem" -down-> "AdditiveGroupElem"
"AdditiveGroupElem" -down-> "ModuleElem{T}"
"ModuleElem{T}" -down-> "FPModuleElem{T}"
"ModuleElem{T}" -down-> "MatElem{T}"

' Non-commutative rings
"SetElem" -down-> "NCRingElem"
"NCRingElem" -down--> "MatAlgElem{T}"
"NCRingElem" -down--> "NCPolyElem{T}"
"NCRingElem" -down-> "RingElem"

' Rings
"RingElem" -down-> "PolyElem{T}"
"PolyElem{T}" -down-> "LaurentPolyElem{T}"
"PolyElem{T}" -down-> "SeriesElem{T}"
"RingElem" -down-> "ResElem{T}"
"RingElem" -down---> "FieldElem"

' Fields
"FieldElem" -down-> "ResFieldElem{T}"
"FieldElem" -down-> "FracElem{T}"
"FieldElem" -down-> "NumFieldElem{T}"
"NumFieldElem{T}" -down-> "SimpleNumFieldElem{T}"
"FieldElem" -down-> "FinFieldElem"

hide members
hide circle

@enduml
"""
```
```@example elements
elements #hide
```
