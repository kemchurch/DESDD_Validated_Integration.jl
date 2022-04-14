function main(filename)
    D = load(filename)
    ppi = D["ppi"]
    δ₋ = D["δ₋"]
    T = D["T"]
    n₊ = D["n₊"]
    para = D["para"]
    r_∞ = D["r_∞"]
    ipara = D["ipara"]
    r_star = D["r_star"]
    rounding_YZ = D["rounding_YZ"]
    r_0 = D["r_0"]
    𝐤 = D["𝐤"]
    k = D["k"]
    m = D["m"]
    u_jℓ = D["u_jℓ"]
    δ = D["δ"]
    Φ₍₀₎ = D["Φ₍₀₎"]
    kΦ = D["kΦ"]
    kΦ_low = D["kΦ_low"]
    Φ,Φ_function,δ,r_δ,r_C⁰ = _main( u_jℓ, δ, ppi, T, Φ₍₀₎, δ₋, para, ipara, n₊, r_star, r_0, r_∞, 𝐤, k, m, kΦ, kΦ_low, rounding_YZ )
    return Φ,Φ_function,δ,r_δ,r_C⁰
end

function _main(u_jℓ, δ,ppi, T, Φ₍₀₎, δ₋, para, ipara, n₊, r_star, r_0, r_∞, 𝐤, k, m, kΦ, kΦ_low, rounding_YZ)
    setprecision(IntervalArithmetic.Interval,512)
    setprecision(BigFloat,512)
    # Prove solutions -----------------------------------------------------
    # Initialize containers.
    δ_val = Interval{T}[]
    U_jℓ = Sequence{CartesianPower{Chebyshev}, Vector{T}}[]
    Φ = Sequence{CartesianPower{Chebyshev}, Vector{Interval{T}}}[]
    Φ_low_accuracy = Sequence{CartesianPower{Chebyshev}, Vector{Interval{T}}}[]
    Φ_function = Array{Any}(undef,n₊)   # Yes, this is inefficient.
    r_C⁰ = Interval{T}[]
    r_δ = Interval{T}[]
    # First step 
    # # -- preparation.
    printstyled("::: Starting proof 1 :::\n",color = :yellow)
    push!(Φ,interval.(Φ₍₀₎))
    push!(Φ_low_accuracy,interval.(Φ₍₀₎))
    Φ_function[1] = convert_Chebyshev_to_callable_function(Φ[1])
    s_jℓ = make_s_grid(k[1],m[1],ppi;T);
    push!(U_jℓ,convert_matrix_to_interp(u_jℓ[1],mid.(ppi)))
    println(": Entering Newton :")
    u_jℓ[1],U_jℓ[1],δ[1],res = Newton((u_jℓ[1],U_jℓ[1]),δ[1],mid.(s_jℓ),mid.(Φ[1]),para,mid.(ppi),1E-20,12;δ₋)
    # # -- conversion to intervals, evaluation of zero-finding problem and approximate derivative.
    iu = (Interval.(u_jℓ[1]),Interval.(U_jℓ[1]));    iδ = Interval(δ[1]);   iΦ = Interval.(Φ[1]);   iδ₋ = interval(δ₋);
    iu, iδ =  modify_candidate_zero_rightside_only(iu,iδ,para;desired_tolerance=res)
    println(": Computing G, DG, A :")
    _,G = Gₚ(iu[1],iu[2],s_jℓ,iδ,component(iΦ,1),ipara,ppi;δ₋=iδ₋);
    DG = DGₚ(iu,iδ,s_jℓ,component(iΦ,1:2),ipara,ppi;δ₋=iδ₋);
    A = Interval.(compute_A(DG,Float64));
    # # -- Evaluate radii polynomial.
    println(": Starting evaluation of the bounds, radii polynomials. :")
    _, _, _, _, _, _, _, _, _, _, C0_err, δ_err = radpol(G,A,DG,iu,iδ,iΦ,ipara,s_jℓ,ppi,r_star[1],r_0[1],r_∞[1];δ₋=iδ₋,rounding_YZ);
    push!(r_C⁰,C0_err)
    push!(r_δ,δ_err)
    push!(δ_val, iδ + interval(-1,1)*r_δ[1])
    # # -- Interpolate, get ready for step 2.
    println(": Starting high-accuracy interpolation. :")
    _,_,Φ_interp_high_accuracy = interpolate_solution_tight_D012(interval.(U_jℓ[1]),interval.(δ[1]),s_jℓ,Φ[1],r_C⁰[1],r_δ[1],kΦ[1],interval(kΦ[1]),𝐤,ipara,ppi;δ₋=iδ₋,sdiv=10,check_large_coeffs=0)
    println(": Starting low-accuracy interpolation. :")
    _,_,Φ_interp_low_accuracy = interpolate_solution_tight_D012(interval.(U_jℓ[1]),interval.(δ[1]),s_jℓ,Φ[1],r_C⁰[1],r_δ[1],kΦ_low[1],interval(kΦ_low[1]),𝐤,ipara,ppi;δ₋=iδ₋,sdiv=10,max_N=kΦ_low[1],check_large_coeffs=0)
    push!(Φ,copy(Φ_interp_high_accuracy))
    push!(Φ_low_accuracy,copy(Φ_interp_low_accuracy))
    println(": Building hybrid enclosure. :")
    Φ_function[2] = t -> evaluation_hybrid_enclosure(t,convert_Chebyshev_to_callable_function(component(Φ[2],1)),
                                                        δ_val[1],Φ[1],ipara,all_derivatives(u_jℓ[1],ipara,Φ[1],δ_val[1];δ₋=copy(iδ₋));n_derivatives=𝐤,δ₋=copy(iδ₋))
    println(": Verifying monotonicity of time lag. :")
    check_lag_monotonicity(ipara,Φ[1],iu[2],s_jℓ,δ_val[1],r_C⁰[1])
    println(": Calculations for proof 1 complete. :")
    # Continue....
    for n=2:n₊
        printstyled("::: Starting proof $n :::\n",color = :yellow)
        s_jℓ_n = make_s_grid(k[n],m[n],ppi;T)
        push!(U_jℓ,convert_matrix_to_interp(u_jℓ[n],mid.(ppi)))
        local C0_err, δ_err, Φ_interp_high_accuracy, Φ_interp_low_accuracy
        println(": Entering Newton with low-accuracy interpolant (faster) :")
        u_jℓ[n],U_jℓ[n],δ[n],_ = Newton((u_jℓ[n],U_jℓ[n]),δ[n],mid.(s_jℓ_n),mid.(Φ_low_accuracy[n]),para,mid.(ppi),1E-30,12;δ₋=δ[n-1]) 
        Newton_tol = min(sup(r_C⁰[n-1]),sup(r_δ[n-1]))/(10*m[n])                                                    
        println(": Entering Newton with mixed accuracy interpolant (G - high accuracy / DG - low accuracy) :")                        
        u_jℓ[n],U_jℓ[n],δ[n] = Newton_multiple_Phi((u_jℓ[n],U_jℓ[n]),δ[n],mid.(s_jℓ_n),mid.(Φ[n]),mid.(Φ_low_accuracy[n]),para,mid.(ppi),Newton_tol,12;δ₋=δ[n-1])                                     
        iu_n = (Interval.(u_jℓ[n]),Interval.(U_jℓ[n]));    iδ_n = Interval(δ[n]);
        iu_n,iδ_n = modify_candidate_zero(iu_n,iδ_n,ipara,δ_val[n-1]);
        println(": Computing G, DG, A :")
        _,G_n = Gₚ(iu_n[1],iu_n[2],s_jℓ_n,iδ_n,component(Φ[n],1),ipara,ppi;δ₋=δ_val[n-1],check_edges=1,rigorous_error_control=1);
        DG_n = DGₚ(iu_n,iδ_n,s_jℓ_n,component(Φ_low_accuracy[n],1:2),ipara,ppi;δ₋=δ_val[n-1],rigorous_error_control=1);
        A_n = Interval.(compute_A(DG_n,Float64));
        if n<n₊
            C0_err, δ_err, Φ_interp_high_accuracy, Φ_interp_low_accuracy, Φ_function[n+1] = proof_section_radpol_interp_monotonicity( G_n,A_n,DG_n,iu_n,iδ_n,Φ[n],Φ_function[n],ipara,s_jℓ_n,ppi,r_star[n],r_0[n],r_∞[n],δ_val[n-1],kΦ[n],kΦ_low[n],𝐤,rounding_YZ )
        else
            C0_err, δ_err, Φ_interp_high_accuracy, Φ_interp_low_accuracy = proof_section_radpol_interp_monotonicity_nohybrid( G_n,A_n,DG_n,iu_n,iδ_n,Φ[n],Φ_function[n],ipara,s_jℓ_n,ppi,r_star[n],r_0[n],r_∞[n],δ_val[n-1],kΦ[n],kΦ_low[n],𝐤,rounding_YZ )
        end
        push!(δ_val, iδ_n + interval(-1,1)*δ_err)
        push!(r_C⁰,C0_err)
        push!(r_δ,δ_err)
        push!(Φ,Φ_interp_high_accuracy)
        push!(Φ_low_accuracy,Φ_interp_low_accuracy)
        println(": Calculations for proof $n complete. :")
    end
    return Φ,Φ_function,δ_val,r_δ,r_C⁰
end