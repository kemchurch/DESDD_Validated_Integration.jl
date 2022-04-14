# DESDD_Validated_Integration
DE_SDD_Validated_Integration.jl is a Julia package used to complete the computer-assisted proofs in Section 5 of the paper Validated integration of differential equations with state-dependent delay. 

The package exports a single function called `main`, which returns 
- `\Phi<TAB>` : An array of Chebyshev sequences corresponding to the initial condition (`\Phi[1]`) and several of its derivatives, and the the solutions of four implicit step (`\Phi[2]` to `\Phi[5]`), rescaled to the domain [-1,1], with validation radius included in zeroth coefficients.
- `\Phi<TAB>_function` : An array of functions, such that `\Phi_function[n](t)`for n=1,...,4 returns the value of the hybrid enclosure of `\Phi[n]` and several of its derivatives at argument `t::Interval{Real}` and several of its derivatives.
- `\delta<TAB>` : Vector of validated interval step sizes.
- `r_\delta<TAB>` : Vector of validation radii for the step sizes.
- `r_C\^0<TAB>` : Vector of validation radii for the solution of the DE-SDD.

Each file in the folder [data](https://github.com/kemchurch/DESDD_Validated_Integration.jl/tree/main/data)
