using Base.Test, Nemo

# write your own tests here

pkgdir = Pkg.dir("Nemo")
pwd = "$pkgdir/src"
push!(DL_LOAD_PATH, "$pwd/../src/lib")

Nemo.Rings.Test.test_all()

Nemo.Fields.Test.test_all()

