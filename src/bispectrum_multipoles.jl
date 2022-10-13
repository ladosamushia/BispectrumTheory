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
    return hcubature(B5D, [0, 0], [1, 2*pi])[1]
end