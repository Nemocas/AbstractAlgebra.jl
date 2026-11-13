###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::UniversalRing{MPoly{EuclideanRingResidueRingElem{BigInt}}})
   return rand(R, 0:4, 0:10, -10:10)
end

function ConformanceTests.generate_element(R::UniversalRing{<:LaurentMPolyWrap{EuclideanRingResidueRingElem{BigInt}}})
   return rand(R, 0:4, -10:10, -10:10)
end
