"""
    This file will contain functions for computing the tree level Bispectrum
    model.
"""

using Interpolations

"""
    solve_triangular_geometry(μ1, ϕ, k1, k2, k3)

    Find all relevant angles of the triangle. μ_i is the cosine of the angle 
    of the ki vector with respect to the z-axis. μij is the cosine of the
    angle between vectors ki and kj.

    Input:
    -μ1, ϕ, k1, k2, k3
    Output:
    -μ2, μ3, μ12, μ23, μ31
"""
function solve_triangular_geometry(μ1, ϕ, k1, k2, k3)
    μ12 = (k1^2 + k2^2 - k3^2)/(2*k1*k2)
    μ31 = (k3^2 + k1^2 - k2^2)/(2*k3*k1)
    μ23 = (k2^2 + k3^2 - k1^2)/(2*k2*k3)
    μ3 = - μ1*μ12 - sqrt(1 - μ1^2)*sqrt(1 - μ12^2)*cos(ϕ)
    μ2 = -(μ1*k1 + μ3*k3)/k2
    return μ2, μ3, μ12, μ23, μ31
end

"""
    tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f)

    compute tree level bispectrum in 3D

    -Input:
    -μ1, ϕ::Float64: angles describing bispectrum orientation (see. Gagrani &
     Samushia). ϕ in radians.
    - k1, k2, k3::Array{Float64,1}: wavenumbers describing the triangle shape
    - b1, b2, f::Float64: bias parameters and the growth rete. 
    - pk::Interpolations.GriddedInterpolation: interpolator for the power
    spectrum
    -Return:
    -Bisp::Array{Float64,1}: bispectrum at those triangles. Same size as k1,k2
    k3
"""
function tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)

    pk1 = pk(k1)
    pk2 = pk(k2)
    pk3 = pk(k3)

    μ2, μ3, μ12, μ23, μ31 = solve_triangular_geometry(μ1, ϕ, k1, k2, k3)
 
    Z1k1 = b1 + f*μ1^2
    Z1k2 = b1 + f*μ2^2
    Z1k3 = b1 + f*μ3^2
 
    F12 = 5/7 + μ12/2*(k1/k2 + k2/k1) + 2/7*μ12^2
    F23 = 5/7 + μ23/2*(k2/k3 + k3/k2) + 2/7*μ23^2
    F31 = 5/7 + μ31/2*(k3/k1 + k1/k3) + 2/7*μ31^2
 
    G12 = 3/7 + μ12/2*(k1/k2 + k2/k1) + 4/7*μ12^2
    G23 = 3/7 + μ23/2*(k2/k3 + k3/k2) + 4/7*μ23^2
    G31 = 3/7 + μ31/2*(k3/k1 + k1/k3) + 4/7*μ31^2
 
    Z2k12 = b2/2. + b1*F12 + f*μ3^2*G12
    Z2k12 -= f*μ3*k3/2*(μ1/k1*Z1k2 + μ2/k2*Z1k1)
    Z2k23 = b2/2. + b1*F23 + f*μ1^2*G23
    Z2k23 -= f*μ1*k1/2*(μ2/k2*Z1k3 + μ3/k3*Z1k2)
    Z2k31 = b2/2. + b1*F31 + f*μ2^2*G31
    Z2k31 -= f*μ2*k2/2*(μ3/k3*Z1k1 + μ1/k1*Z1k3)

    Bi = Z2k12*Z1k1*Z1k2*pk1*pk2
    Bi += Z2k23*Z1k2*Z1k3*pk2*pk3
    Bi += Z2k31*Z1k3*Z1k1*pk3*pk1
 
    return 2*Bi 
end