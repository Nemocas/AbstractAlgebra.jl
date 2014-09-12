using Base.Test, Nemo

# write your own tests here

pwd = chomp(readall(`pwd`))
push!(DL_LOAD_PATH, "$pwd/../src/lib")

Nemo.Rings.Test.test_all()

Nemo.Fields.Test.test_all()

