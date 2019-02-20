include("generic/FreeModule-test.jl")
include("generic/ModuleMorphism-test.jl")
include("generic/Submodule-test.jl")

function test_modules()
   test_free_module()
   test_module_morphism()
   test_submodule()
end
