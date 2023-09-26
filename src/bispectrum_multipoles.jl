"""
    This file will contain functions that integrate 5D bispectrum to get
    multipoles.
"""

using Cubature

include("./TreeLevelBispectrum.jl")

"""
    B00(k1, k2, k3, b1, b2, f, pk, bk)

    Input:
    -k1, k2, k3::Float64 - wavenumbers for the bispectrum shape
    -b1, b2, f::Float64 - bias parameters
    -pk::interpolator - interpolator for the pk
    -bk::function - function to compute bispectrum
"""
function B00(k1, k2, k3, b1, b2, f, pk, bk)
    B5D = x -> bk(x[1], x[2], k1, k2, k3, b1, b2, f, pk)
    return hcubature(B5D, [-1, 0], [1, 2*pi], reltol=1e-8)[1]/8/π
end

function B00_itp(k, b1, b2, f, pk, bk)
    N = length(k)
    B00_grid = zeros(N, N, N)
    for i1 in 1:N, i2 in i1:N, i3 in i2:N
        Bth = B00(k[i1], k[i2], k[i3], b1, b2, f, pk, bk)
        B00_grid[i1,i2,i3] = Bth
        B00_grid[i1,i3,i2] = Bth
        B00_grid[i2,i1,i3] = Bth
        B00_grid[i2,i3,i1] = Bth
        B00_grid[i3,i1,i2] = Bth
        B00_grid[i3,i2,i1] = Bth
    end
    print(B00_grid)
    itp = interpolate((k1, k2, k3), B00_grid, BSpline(Linear()))
    return itp
end
