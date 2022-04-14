module DESDD_Validated_Integration

using JLD2, RadiiPolynomial

include("main_functions.jl")
include("core_proof_functions.jl")
include("chebyshev_helper.jl")
include("vector_fields.jl")
include("piecewise_interpolation.jl")

export main

end