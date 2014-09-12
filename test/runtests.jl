using Base.Test, Nemo

# write your own tests here

pwd = chomp(readall(`pwd`))
println(pwd)
push!(DL_LOAD_PATH, "$pwd/lib")

Nemo.Rings.Test.test_all()

Nemo.Fields.Test.test_all()

