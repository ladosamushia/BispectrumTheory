"""
    Testing functions in the bispectrum_multipoles.jl file
"""

using Test

include("../src/bispectrum_multipoles.jl")

@testset "bispectrum_multipoles" begin
    pk = interpolate(([0, 1, 2, 3],), [1, 1, 1, 1], Gridded(Linear()))
    bk = (μ1, ϕ, k1, k2, k3, b1, b2, f, pk) -> μ1^2*ϕ + μ1*cos(ϕ)
    @test isapprox(B00(1, 1, 1, 1, 1, 1, pk, bk), (2*π)^2/3/2)
    bk = (μ1, ϕ, k1, k2, k3, b1, b2, f, pk) -> μ1^2*cos(ϕ/4)
    @test isapprox(B00(1, 1, 1, 1, 1, 1, pk, bk), 4/3)
end;