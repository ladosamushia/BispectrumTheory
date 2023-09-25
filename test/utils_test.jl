using Test

include("../src/utils.jl")

@testset "solve_triangular_geometry" begin
    computation = solve_triangular_geometry(0,0,1,1,1)
    expectation = [sqrt(3)/2, -sqrt(3)/2, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
    computation = solve_triangular_geometry(1,0,1,1,1)
    expectation = [-1/2, -1/2, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
    computation = solve_triangular_geometry(sqrt(3)/2,0,1,1,1) 
    expectation = [0, -sqrt(3)/2, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
    computation = solve_triangular_geometry(1/2,0,1,1,1)
    expectation = [1/2, -1, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
    computation = solve_triangular_geometry(1,pi/3,1,1,1)
    expectation = [-1/2, -1/2, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
    computation = solve_triangular_geometry(0,pi/2,1,1,1)
    expectation = [0, 0, 1/2, 1/2, 1/2]
    @test all(abs.(computation .- expectation) .< 1e-6)
end

@testset "kkk_grid" begin
    kmin = 0.05; kstep = 0.1; Nk = 3
    computation = kkk_grid(kmin, kstep, Nk)
    expectation = [0.05 0.15 0.15 0.25 0.25 0.25 0.25; 0.05 0.15 0.15 0.15 0.25 0.25 0.25; 0.05 0.05 0.15 0.15 0.05 0.15 0.25]
    @test all(abs.(computation .- expectation) .< 1e-6)
end