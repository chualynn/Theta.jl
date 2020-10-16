import Theta: parity_char, remainder_char

@testset "Characteristics" begin
    @test parity_char([[0,0],[0,1]]) == 0
    @test parity_char([[1,1],[0,1]]) == 1
    @test remainder_char([[0,1,2],[3,1,5]]) == [[0,1,0],[1,1,1]]
    @test Set(theta_char(1)) == Set([[[0],[0]], [[0],[1]], [[1],[0]], [[1],[1]]])
    @test Set(even_theta_char(1)) == Set([[[0],[0]], [[0],[1]], [[1],[0]]])
    @test odd_theta_char(1) == [[[1],[1]]]
    @test check_azygetic([[[1,0,1,0], [1,0,1,0]], [[0,0,0,1], [1,0,0,0]], [[0,0,1,1], [1,0,1,1]]]) == true
end
