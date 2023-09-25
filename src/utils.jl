"""
    Utility functions used elsewhere in the code.
"""

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
    μ3 = - μ1*μ31 - sqrt(1 - μ1^2)*sqrt(1 - μ31^2)*cos(ϕ)
    μ2 = -(μ1*k1 + μ3*k3)/k2
    return μ2, μ3, μ12, μ23, μ31
end


"""
    kkk_grid(kmin, kstep, Nk)

    Construct a 3D grid of k1, k2, k3 so that k1 > k2 > k3 and the triangular
    condition is satisfied.

    Input:
    -kmin::Float64: minimum value of k1, k2, k3
    -kstep::Float64: size of the k bins
    -Nk::Int: number of k bins

    Output:
    -kkk::Array{Float64,2}: 3xNk numbers with bin centers for k1, k2, k3
    triplets.
"""
function kkk_grid(kmin, kstep, Nk)
    kgrid = collect(range(kmin, step=kstep, length=Nk))
    Ntri = 0
    for k1 in kgrid, k2 in kgrid, k3 in kgrid
        if k1 <= k2 <= k3 && k1 + k2 >= k3
            Ntri += 1
        end
    end
    kkk = zeros(3, Ntri)
    counter = 1
    for k1 in kgrid, k2 in kgrid, k3 in kgrid
        if k1 <= k2 <= k3 && k1 + k2 >= k3
            kkk[:,counter] = [k1; k2; k3]
            counter += 1
        end
    end
    return kkk
end