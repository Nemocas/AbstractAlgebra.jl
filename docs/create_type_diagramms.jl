using Kroki

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
"Ring" -down---> "Field"
"Ring" -down-> "PolyRing{T}"
"PolyRing{T}" -down-> "LaurentPolynomialRing{T}"
"PolyRing{T}" -down-> "SeriesRing{T}"
"Ring" -down-> "ResRing{T}"

' Fields
"Field" -down-> "ResField{T}"
"NumField{T}" -down-> "SimpleNumField{T}"
"Field" -down-> "FracField{T}"
"Field" -down-> "NumField{T}"
"Field" -down-> "FinField"

hide members
hide circle

@enduml
"""

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
"PolyElem{T}" -down-> "SeriesElem{T}"
"PolyElem{T}" -down-> "LaurentPolyElem{T}"
"RingElem" -down-> "ResElem{T}"
"RingElem" -down---> "FieldElem"

' Fields
"FieldElem" -down-> "FinFieldElem"
"FieldElem" -down-> "ResFieldElem{T}"
"FieldElem" -down-> "NumFieldElem{T}"
"NumFieldElem{T}" -down-> "SimpleNumFieldElem{T}"
"FieldElem" -down-> "FracElem{T}"

hide members
hide circle

@enduml
"""


open(joinpath(@__DIR__, "src", "parents_diagram.svg"), "w") do io
  write(io, sprint(show, "image/svg+xml", parents))
end

open(joinpath(@__DIR__, "src", "elements_diagram.svg"), "w") do io
  write(io, sprint(show, "image/svg+xml", elements))
end
