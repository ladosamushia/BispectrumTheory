"""
    This file will contain tests associated with the TreeLevelBispectrum.jl
    file.
"""

using Test

include("../src/TreeLevelBispectrum.jl")

# Do some limiting values that are easy to compute by hand
#@testset "TreeLevelBispectrum" begin
    pk = interpolate(([0, 1, 2, 3],), [1, 1, 1, 1], Gridded(Linear()))

    #μ1 = 0; ϕ = 0; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 1; b2 = 1; f = 0
    #bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    #@test isapprox(bk, 3*2*25/14)

    μ1 = 0; ϕ = 0; k1 = 1.0; k2 = 1.0; k3 = 1.0; b1 = 0; b2 = 0; f = 1
    bk = tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)
    #@test isapprox(bk, 1/28)
#end;