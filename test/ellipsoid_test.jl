import Theta: ellipsoid_pointwise, ellipsoid_uniform

@testset "Ellipsoid" begin
    T = [1 0; 0 1];
    shift = [0; 0];
    @test Set(ellipsoid_pointwise(T, 2, shift)) == Set([[0,-1], [-1,0], [1,0], [0,1], [0,0]])
    @test issubset(ellipsoid_pointwise(T, 2, shift), ellipsoid_uniform(T, 2))
end
