"""
    This file will contain tests associated with the TreeLevelBispectrum.jl
    file.
"""

using Test

include("../src/TreeLevelBispectrum.jl")

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

# Do some limiting values that are easy to compute by hand
@testset "TreeLevelBispectrum" begin
    pk = interpolate(([0, 1, 2, 3],), [1, 1, 1, 1], Gridded(Linear()))

    μ1 = 0; ϕ = 0; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 1; b2 = 1; f = 0
    bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    @test isapprox(bk, 3*2*25/14)

    μ1 = 0; ϕ = 0; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 1; b2 = 0; f = 1
    bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    @test isapprox(bk, 201/8)

    μ1 = 1/2; ϕ = 0; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 1; b2 = 0; f = 1
    bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    @test isapprox(bk, 11415/448)

    μ1 = 1; ϕ = pi/2; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 1; b2 = 0; f = 1
    bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    @test isapprox(bk, 11415/448)
end;