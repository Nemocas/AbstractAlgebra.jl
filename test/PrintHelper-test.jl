@testset "pluralize" begin
  @test pluralize("basis") == "bases"
  @test pluralize("child") == "children"
  @test pluralize("generator") == "generators"
  @test pluralize("indeterminate") == "indeterminates"
  @test pluralize("matrix") == "matrices"
  @test pluralize("maximum") == "maxima"
  @test pluralize("polyhedron") == "polyhedra"
  @test pluralize("tetrahedron") == "tetrahedra"
  @test pluralize("variable") == "variables"
  @test pluralize("variety") == "varieties"
  @test pluralize("vertex") == "vertices"
end

@testset "ItemQuantity" begin
  @test string(ItemQuantity(0, "generator")) == "0 generators"
  @test string(ItemQuantity(1, "generator")) == "1 generator"
  @test string(ItemQuantity(2, "generator")) == "2 generators"
  @test string(ItemQuantity(0, "ox", "oxen")) == "0 oxen"
  @test string(ItemQuantity(1, "ox", "oxen")) == "1 ox"
  @test string(ItemQuantity(2, "ox", "oxen")) == "2 oxen"
end
