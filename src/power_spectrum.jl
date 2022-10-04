"""
    Codes related to the power spectrum of matter and galaxies.
"""

using DelimitedFiles
using Interpolations

"""
    init_linear_pk(filename)

    This function reads linear matter power spectrum from a file and returns an
    interpolator.

    Input:
    - filename::String - name of the file with the linear power spectrum
    
    Output:
    - linPK:: - interpolation object for the linear power spectrum
"""

function init_linear_pk(filename)
    kPk = readdlm(filename, comment="#")
