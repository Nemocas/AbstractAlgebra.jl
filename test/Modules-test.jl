function rand_module(R::AbstractAlgebra.Ring, vals...)
   rk = rand(0:5)
   M = FreeModule(R, rk)
   levels = rand(0:3)
   for i = 1:levels
      if ngens(M) == 0
         break
      end
      G = [rand(M, vals...) for i in 1:rand(1:ngens(M))]
      S, f = sub(M, G)
      if rand(1:2) == 1
         M, f = quo(M, S)
      else
         M = S
      end
   end
   return M
end


include("generic/FreeModule-test.jl")
include("generic/ModuleHomomorphism-test.jl")
include("generic/Submodule-test.jl")
include("generic/QuotientModule-test.jl")
include("generic/DirectSum-test.jl")
include("generic/Module-test.jl")
include("generic/InvariantFactorDecomposition-test.jl")
