@testset "error printing" begin
    b = IOBuffer()
    Base.showerror(b, NotImplementedError(:foo, 1, 2))
    Base.showerror(b, NotInvertibleError(ZZ(0)))
    @test length(take!(b)) > 1
end
