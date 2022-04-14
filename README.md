# DESDD_Validated_Integration.jl
DE_SDD_Validated_Integration.jl is a Julia package used to complete the computer-assisted proofs in Section 5 of the paper Validated integration of differential equations with state-dependent delay. 

The package exports a single function called `main`, which returns 
- `Φ` : An array of Chebyshev sequences corresponding to the initial condition (`Φ[1]`) and several of its derivatives, and the the solutions of four implicit step (`Φ[2]` to `Φ[5]`), rescaled to the domain [-1,1], with validation radius propagated forward by inclusion in the order zero term.
- `Φ_function` : An array of functions, such that `Φ_function[n](t)`for n=1,...,4 returns the value of the hybrid enclosure of `Φ[n+1]` and several of its derivatives at argument `t::Interval{T} where T<:Real`.
- `δ` : Vector of validated interval step sizes.
- `r_δ` : Vector of validation radii for the step sizes.
- `r_C⁰` : Vector of C⁰ enclosures for the solution of the DE-SDD.

Each file in the folder [data](https://github.com/kemchurch/DESDD_Validated_Integration.jl/tree/main/data) has a name formatted in the style `P_X_initialcondition_Y.jld2`. The symbol `X` references the parameter set, while Y references the initial condition (see Section 5.2 of the paper). These files contain variables which include numerically computed solutions of the relevant initial-value problems, as well as parameters required for the computer-assisted proofs.

## Installation
The package requires Julia v1.6 or higher to be installed. Clone this repository, and activate/instantiate the package as follows:
```julia
import Pkg
Pkg.activate("path/to/package") # edit the path accordingly
Pkg.instantiate()
```

## Usage
To complete the proofs described in Section 5 of the paper, execute the following for each file in the data folder. Be aware that the proofs can require a lot of time to complete.

```julia
using DESDD_Validated_Integration
filename = "path/to/data/file.jld2" # edit the path accordingly (cf. folder data)
Φ,Φ_function,δ,r_δ,r_C⁰ = main(filename)
```
