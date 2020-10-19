@testset "Riemann matrix" begin
    τ = random_siegel(2);
    R = RiemannMatrix(τ);
    @test R.τ == τ
    @test R.g == 2
    @test R.X == real(τ)
    @test R.Y == imag(τ)
    @test size(R.ellipsoid)[1] == 5
    @test R.T[2,1] == 0
end
