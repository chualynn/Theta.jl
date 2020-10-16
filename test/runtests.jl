using Test
using Theta
using LinearAlgebra

@testset "Theta" begin
    include("theta_test.jl")
    include("lattice_test.jl")
    include("characteristics_test.jl")
    include("accola_test.jl")
end
