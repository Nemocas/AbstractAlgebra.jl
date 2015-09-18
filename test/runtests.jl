using Nemo, Base.Test

pkgdir = Pkg.dir("Nemo")
pwd = "$pkgdir/src"

push!(Libdl.DL_LOAD_PATH, "$pwd/../src/lib")

Nemo.Test.test_all()
