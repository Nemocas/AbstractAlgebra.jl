function test_submodule()
   M = FreeModule(QQ, 2)

   for iter = 1:200
      M = FreeModule(QQ, 5)
      nsubs = rand(0:5)
      subs = [Submodule(M, [rand(M, -10:10)]) for i in 1:nsubs]
      N = Submodule(M, subs)
   end
end

