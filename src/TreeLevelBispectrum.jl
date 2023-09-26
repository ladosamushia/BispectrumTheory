"""
    This file will contain functions for computing the tree level Bispectrum
    model.
"""

using Interpolations

include("./utils.jl")

"""
    tree_level_bk(μ1, ϕ, k1, k2, k3, b1, b2, f, pk)

    compute tree level bispectrum in 3D

    -Input:
    -μ1, ϕ::Float64: angles describing bispectrum orientation (see. Gagrani &
     Samushia). ϕ in radians.
    - k1, k2, k3::Array{Float64,1}: wavenumbers describing the triangle shape
    - b1, b2, f::Float64: bias parameters and the growth rete. 
    - pk::Interpolations.GriddedInterpolation: interpolator for the power
    spectrum
    -Return:
    -Bi::Array{Float64,1}: bispectrum at those triangles. Same size as k1,k2
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

"""
    tree_level_b00(k1, k2, k3, b1, b2, f, pk)

    compute tree level bispectrum monopole analitically

    -Input:
    - k1, k2, k3::Array{Float64,1}: wavenumbers describing the triangle shape
    - b1, b2, f::Float64: bias parameters and the growth rete. 
    - pk::Interpolations.GriddedInterpolation: interpolator for the power
    spectrum
    -Return:
    -Bi::Array{Float64,1}: bispectrum at those triangles. Same size as k1,k2
    k3
"""
function tree_level_b00(k1, k2, k3, b1, b2, f, pk)
    
    pk1 = pk(k1)
    pk2 = pk(k2)
    pk3 = pk(k3)

    _, _, μ12, μ23, μ31 = solve_triangular_geometry(0, 0, k1, k2, k3)
    ν12 = sqrt(1 - μ12^2)
    ν23 = sqrt(1 - μ23^2)
    ν31 = sqrt(1 - μ31^2)
 
    F12 = 5/7 + μ12/2*(k1/k2 + k2/k1) + 2/7*μ12^2
    F23 = 5/7 + μ23/2*(k2/k3 + k3/k2) + 2/7*μ23^2
    F31 = 5/7 + μ31/2*(k3/k1 + k1/k3) + 2/7*μ31^2
 
    G12 = 3/7 + μ12/2*(k1/k2 + k2/k1) + 4/7*μ12^2
    G23 = 3/7 + μ23/2*(k2/k3 + k3/k2) + 4/7*μ23^2
    G31 = 3/7 + μ31/2*(k3/k1 + k1/k3) + 4/7*μ31^2

    # The following long lines were computed in a WolframMathematica noteboook
    # The notebook is attached to this repo
    Bi = (1/(630*k1*k2*k3))*(k3*pk1*pk2*(21*b2*f^2*k1*k2*(3*μ12^2 + 
    ν12^2) + 
         f^4*k3*(35*μ12^3*(k1 - k2*μ12)*μ31 + 
            15*μ12*(k1 - 2*k2*μ12)*μ31*ν12^2 - 
            3*k2*μ31*ν12^4 + 
            5*μ12^2*(-3*k1 + 4*k2*μ12)*ν12*ν31 - 
            3*(k1 - 4*k2*μ12)*ν12^3*ν31) + 
         105*b1^3*(6*F12*k1*k2 - 
            f*k3*(k2*μ31 - k1*μ12*μ31 + 
               k1*ν12*ν31)) + 
         3*b1*f*(35*b2*k1*k2*(1 + μ12^2 + ν12^2) + 
            f*(14*F12*k1*k2*(3*μ12^2 + ν12^2) + 
               3*f*k3*(-k2*μ31*(2 + μ12^2 + ν12^2)*(5*
    μ12^2 + ν12^2) + 
                  k1*μ12*μ31*(5 + 10*μ12^2 + 6*ν12^2) + 
                  4*k2*μ12*ν12*(1 + μ12^2 + ν12^2)*
    ν31 - k1*ν12*(1 + 6*μ12^2 + 2*ν12^2)*ν31))) + 
         21*b1^2*(15*b2*k1*k2 + 
            f*(10*F12*k1*k2*(1 + μ12^2 + ν12^2) + 
               f*k3*(3*k1*μ12*μ31*(2 + μ12^2 + ν12^2) - 
                  k2*μ31*(3 + 6*μ12^2 + 2*ν12^2) + 
                  4*k2*μ12*ν12*ν31 - 
                  k1*ν12*(2 + 3*μ12^2 + 3*ν12^2)*ν31)))) +
       k2*pk1*pk3*(21*b2*f^2*k1*k3*(3*μ31^2 + ν31^2) + 
         f^4*k2*(35*μ12*μ31^3*(k1 - k3*μ31) + 
            5*μ31^2*(-3*k1 + 4*k3*μ31)*ν12*ν31 + 
            15*μ12*μ31*(k1 - 2*k3*μ31)*ν31^2 - 
            3*(k1 - 4*k3*μ31)*ν12*ν31^3 - 
            3*k3*μ12*ν31^4) + 
         105*b1^3*(6*F31*k1*k3 - 
            f*k2*(k3*μ12 - k1*μ12*μ31 + 
               k1*ν12*ν31)) + 
         21*b1^2*(15*b2*k1*k3 + 
            f*(10*F31*k1*k3*(1 + μ31^2 + ν31^2) - 
               f*k2*(-4*k3*μ31*ν12*ν31 - 
                  3*k1*μ12*μ31*(2 + μ31^2 + ν31^2) + 
                  k3*μ12*(3 + 6*μ31^2 + 2*ν31^2) + 
                  k1*ν12*ν31*(2 + 3*μ31^2 + 3*ν31^2)))) + 
         3*b1*f*(35*b2*k1*k3*(1 + μ31^2 + ν31^2) + 
            f*(14*F31*k1*k3*(3*μ31^2 + ν31^2) + 
               3*f*k2*(4*k3*μ31*ν12*ν31*(1 + μ31^2 + 
    ν31^2) - 
                  k3*μ12*(2 + μ31^2 + ν31^2)*(5*μ31^2 + 
    ν31^2) - k1*ν12*ν31*(1 + 6*μ31^2 + 2*ν31^2) + 
                  k1*μ12*μ31*(5 + 10*μ31^2 + 
                     6*ν31^2))))) - 
      k1*pk2*pk3*(-105*b1^3*(6*F23*k2*k3 - 
            f*k1*(k3*μ12 + k2*μ31)) - 
         21*b2*f^2*k2*k3*(-4*μ12*μ31*ν12*ν31 + μ12^2*
    (3*μ31^2 + ν31^2) + ν12^2*(μ31^2 + 3*ν31^2)) - 
         21*b1^2*(15*b2*k2*k3 + 
            10*f*F23*k2*k3*(μ12^2 + μ31^2 + ν12^2 + 
    ν31^2) - 
            f^2*k1*(k3*(3*μ12^3 + 6*μ12*μ31^2 + 
                  3*μ12*ν12^2 - 4*μ31*ν12*ν31 + 
                  2*μ12*ν31^2) + 
               k2*(6*μ12^2*μ31 + 3*μ31^3 + 
                  2*μ31*ν12^2 - 4*μ12*ν12*ν31 + 
                  3*μ31*ν31^2))) + 
         f^4*k1*(k2*(-12*μ12^3*ν12*ν31*(5*μ31^2 + 
    ν31^2) + 3*μ31*ν12^4*(μ31^2 + 5*ν31^2) - 
               4*μ12*ν12^3*ν31*(9*μ31^2 + 5*ν31^2) + 
               6*μ12^2*μ31*ν12^2*(5*μ31^2 + 
                  9*ν31^2) + 
               5*μ12^4*(7*μ31^3 + 3*μ31*ν31^2)) + 
            k3*(-12*μ12^2*μ31*ν12*ν31*(5*μ31^2 + 
                  3*ν31^2) - 
               4*μ31*ν12^3*ν31*(3*μ31^2 + 
                  5*ν31^2) + μ12^3*(35*μ31^4 + 
                  30*μ31^2*ν31^2 + 3*ν31^4) + 
               3*μ12*ν12^2*(5*μ31^4 + 
                  18*μ31^2*ν31^2 + 5*ν31^4))) - 
         3*b1*f*(35*b2*k2*k3*(μ12^2 + μ31^2 + ν12^2 + 
    ν31^2) + 
            f*(14*F23*k2*k3*(-4*μ12*μ31*ν12*ν31 + 
    μ12^2*(3*μ31^2 + ν31^2) + ν12^2*(μ31^2 + 
                     3*ν31^2)) - 
               3*f*k1*(k3*(-12*μ12^2*μ31*ν12*ν31 + 
                     2*μ12^3*(5*μ31^2 + ν31^2) - 
                     4*μ31*ν12*ν31*(μ31^2 + ν12^2 + 
    ν31^2) + μ12*(μ31^2 + ν31^2)*(5*μ31^2 + 
                        6*ν12^2 + ν31^2)) + 
                  k2*(5*μ12^4*μ31 - 4*μ12^3*ν12*ν31 - 
                     4*μ12*ν12*ν31*(3*μ31^2 + ν12^2 + 
    ν31^2) + μ31*ν12^2*(2*μ31^2 + ν12^2 + 
                        6*ν31^2) + 
                     2*μ12^2*μ31*(5*μ31^2 + 
                        3*(ν12^2 + ν31^2))))))))

    return Bi
end