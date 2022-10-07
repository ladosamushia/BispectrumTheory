"""
    This file will contain functions that integrate 5D bispectrum to get
    multipoles.
"""

using Cubature

include("./TreeLevelBispectrum.jl")

function B00(k1, k2, k3, b1, b2, f, pk)
    B5D = x -> TreeLevelBispectrum(x[1], x[2], k1, k2, k3, b1, b2, f, pk)
    return hcubature(B5D, [0, 1], [0, 2*pi])[1]
end