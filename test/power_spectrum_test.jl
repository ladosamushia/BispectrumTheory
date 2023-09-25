"""
    Testing functions in power_spectrum.jl
"""

using Test

include("../src/power_spectrum.jl")

@testset "init_linear_pk" begin
    pk = init_linear_pk("data/test_pk.txt")
    # Check power spectrum at two randomly chosen nodes from that file
    @test pk(0.12712E-03) == 0.71690E+03
    @test pk(0.23923E+00) == 0.13125E+04
    # Check that interpolation is in between the two values
    @test 0.27384E+05 < pk(0.21000E-01) < 0.27618E+05
end;