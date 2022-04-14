@views function piecewise_to_smooth(U_jℓ, s_jℓ, N_new::Integer, Ni_new::T where {T<:Real}, ppi::T where {T<:Real}; U_error = 0, talkative = 0)
    # Interpolation of piecewise ̄u(t). 
    node = cheb_nodes(N_new, Ni_new, ppi)
    node_shift = (node .+ 1) / 2
    s = 2 * s_jℓ .- 1
    evals = Interval.(zeros(eltype(U_jℓ), N_new + 1))
    for j ∈ 1:N_new+1
        k = 1
        while node[j] > s[end, k]
            k += 1
        end
        if node[j] >= s[1, k] && node[j] <= s[end, k]
            evals[j] = component(U_jℓ, k)(node_shift[j] * 2 / (s_jℓ[end, k] - s_jℓ[1, k]) - (s_jℓ[end, k] + s_jℓ[1, k]) / (s_jℓ[end, k] - s_jℓ[1, k]))[0] + interval(-1, 1) * U_error
        elseif ~(node[j] >= s[end, k]) && node[j] > s[1, k]
            if k == size(s_jℓ)[2]
                evals[j] = component(U_jℓ, k)(node_shift[j] * 2 / (s_jℓ[end, k] - s_jℓ[1, k]) - (s_jℓ[end, k] + s_jℓ[1, k]) / (s_jℓ[end, k] - s_jℓ[1, k]))[0] + interval(-1, 1) * U_error
            else
                if talkative==1
                    println("Could not resolve a node. Incrementing N_new ↦ N_new + 1 = ", N_new + 1, ", and trying again.")
                end
                N_new += 1
                Ni_new += 1
                U!, U_coeffs!, evals!, N_new!, Ni_new! = piecewise_to_smooth(U_jℓ, s_jℓ, N_new, Ni_new, ppi; U_error)
                return U!, U_coeffs!, evals!, N_new!, Ni_new!
            end
        end
    end
    if talkative==1
        println("Interpolation at N_new = ", N_new, " is feasible.")
    end
    U_coeffs = cheb_interp(evals, N_new, Ni_new, ppi)
    U = Sequence(Chebyshev(N_new), U_coeffs)
    return U, U_coeffs, evals, N_new, Ni_new
end

@views function compute_V(U_jℓ, s_jℓ, Φ, δ, para, 𝐤; δ₋ = 1, sdiv = 1)
    # Implements V for epsilon = 0.
    𝐤i = Interval(𝐤)
    segment_norms = Interval.(zeros(eltype(U_jℓ), 𝐤, size(s_jℓ)[2]))
    V = Interval.(zeros(eltype(U_jℓ),𝐤))
    i1 = Interval.(one(eltype(U_jℓ)))
    Δ = s_jℓ[end,1]-s_jℓ[1,1];
    for k = 1:size(s_jℓ)[2]
        for j = 1:sdiv
            ival_cheb = union(-i1 + 2 * (j - 1) / (i1 * sdiv), -i1 + 2 * j / (i1 * sdiv))
            ival_t = s_jℓ[1,k] + Δ*((ival_cheb + 1) / 2)
            eval = abs.(𝐅₈(ival_t, component(U_jℓ, k)(ival_cheb)[0], δ, Φ, para; δ₋, components_Φ = 𝐤+1)[1:𝐤]);
            for m=1:𝐤
                segment_norms[m,k] = max(segment_norms[m,k], eval[m])
            end
        end
    end
    for m=1:𝐤
        V[m] = norm(segment_norms[m,:],Inf)/ (2^(m + 1))
    end
    return V
end

function convert_Chebyshev_to_callable_function(Φ)
    Φ_global = (t -> coefficients(Φ(t)))
    return Φ_global
end

function all_derivatives(u_jℓ, para, Φ::Sequence{CartesianPower{Chebyshev}}, δ; δ₋ = 1)
    derivatives_left = [u_jℓ[1, 1]; 𝐅₈(0, u_jℓ[1, 1], δ, Φ, para, 1:8; δ₋)] .* (@interval(2) .^ (-(0:8)))
    derivatives_right = [u_jℓ[end, end]; 𝐅₈(1, u_jℓ[end, end], δ, Φ, para, 1:8; δ₋)] .* (@interval(2) .^ (-(0:8)))
    return derivatives_left, derivatives_right
end

function all_derivatives(u_jℓ, para, Φ, δ; δ₋ = 1)
    𝐤 = length(Φ(interval(0)))-1;
    derivatives_left = [u_jℓ[1, 1]; 𝐅₈(0, u_jℓ[1, 1], δ, Φ, para; δ₋, components_Φ = 𝐤+1)] .* (@interval(2) .^ (-(0:𝐤)))
    derivatives_right = [u_jℓ[end, end]; 𝐅₈(1, u_jℓ[end, end], δ, Φ, para; δ₋,components_Φ = 𝐤+1)] .* (@interval(2) .^ (-(0:𝐤)))
    return derivatives_left, derivatives_right
end

function interpolate_derivatives(U, δ, Φ::Sequence{CartesianPower{Chebyshev}}, para, N_new::Integer, Ni_new::T where {T<:Real}, number_of_derivatives, ppi::T where {T<:Real}; δ₋ = 1, U_error = 0)
    node = cheb_nodes(N_new, Ni_new, ppi)
    node_shift = (node .+ 1) / 2
    mat_interp = Interval.(zeros(eltype(U), number_of_derivatives, N_new + 1))
    cheb_coeffs = Interval.(zeros(eltype(U), N_new + 1, number_of_derivatives))
    for j = 1:N_new+1
        mat_interp[:, j] = 𝐅₈(node_shift[j], U(node[j])[0], δ, Φ, para, 1:number_of_derivatives; δ₋=δ₋, U_error, order = number_of_derivatives, components_Φ = num_components(Φ)) ./ (2 .^ (1:number_of_derivatives))
    end
    mat_interp = mat_interp'
    for k = 1:number_of_derivatives
        cheb_coeffs[:, k] = cheb_interp(mat_interp[:, k], N_new, Ni_new, ppi)
    end
    all_u = Sequence(Chebyshev(N_new)^(number_of_derivatives + 1), [coefficients(U); reshape(cheb_coeffs, (N_new + 1) * number_of_derivatives)])
    return mat_interp, cheb_coeffs, all_u
end

function interpolate_derivatives(U, δ, Φ, para, N_new::Integer, Ni_new::T where {T<:Real}, number_of_derivatives, ppi::T where {T<:Real}; δ₋ = 1, U_error = 0)
    node = cheb_nodes(N_new, Ni_new, ppi)
    node_shift = (node .+ 1) / 2
    mat_interp = Interval.(zeros(eltype(U), number_of_derivatives, N_new + 1))
    cheb_coeffs = Interval.(zeros(eltype(U), N_new + 1, number_of_derivatives))
    for j = 1:N_new+1
        mat_interp[:, j] = 𝐅₈(node_shift[j], U(node[j])[0], δ, Φ, para; δ₋=δ₋, U_error, components_Φ = length(Φ(interval(0))))[1:number_of_derivatives] ./ (2 .^ (1:number_of_derivatives))
    end
    mat_interp = mat_interp'
    for k = 1:number_of_derivatives
        cheb_coeffs[:, k] = cheb_interp(mat_interp[:, k], N_new, Ni_new, ppi)
    end
    all_u = Sequence(Chebyshev(N_new)^(number_of_derivatives + 1), [coefficients(U); reshape(cheb_coeffs, (N_new + 1) * number_of_derivatives)])
    return mat_interp, cheb_coeffs, all_u
end

@views function hybrid_enclosure(t::Interval{T} where T<:Real, U_t, U_t_left, U_t_right, δ, Φ, para, derivatives_at_endpoints; n_derivatives = 8, δ₋=1)
    # Note: error terms can be implicitly passed forward by including them in U_t, U_t_left and U_t_right, and the evaluations of the (bootstrap) endpoint derivatives.
    # Note(2): In all cases, t is in the scaling of [-1,1]. U_t = U((t+1)/2), while U_t_left and U_t_right are U((t₋+1)/2) and U((t₊+1)/2), where t₋ and t₊ are defined
    # such that t = ∪{t₋,t₊} (union disjoint aside from one point) and t₋∩t₊ = t∩{-1,1} is a singleton.
    derivatives_left, derivatives_right = derivatives_at_endpoints
    h_enc_inf = zeros(eltype(U_t), n_derivatives + 1)
    h_enc_sup = zeros(eltype(U_t), n_derivatives + 1)
    inf_t = inf(t)
    sup_t = sup(t)
    for n = 0:n_derivatives
        if inf_t <= -1
            f8 = Interval.(factorial.(0:n_derivatives))
            if sup_t <= -1
                h_enc_inf[n+1] = (sum((derivatives_left[n+1:end] .* ((t + 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
                h_enc_sup[n+1] = (sum((derivatives_left[n+1:end] .* ((t + 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
            else
                t_left = interval(inf_t, -1)
                h_enc_inf[n+1] = (sum((derivatives_left[n+1:end] .* ((t_left + 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
                if sup_t <= 1
                    t_right = interval(-1,sup_t);
                    if n>=1
                        h_enc_sup[n+1] = 𝐅₈((t_right+1)/2, U_t_right, δ, Φ, para,; δ₋)[n] / @interval(2)^n;
                    else
                        h_enc_sup[n+1] = U_t_right;
                    end
                else
                    error("An input (somewhere) is an interval that contains [-1,1]. This could be a result of error propagation. I can't handle that.")
                end
            end
        elseif inf_t <= 1
            if sup_t <= 1
                if n>=1
                    h_enc_inf[n+1] = 𝐅₈((t+1)/2, U_t, δ, Φ, para; δ₋)[n] / @interval(2)^n;
                    h_enc_sup[n+1] = 𝐅₈((t+1)/2, U_t, δ, Φ, para; δ₋)[n] / @interval(2)^n;
                else
                    h_enc_inf[n+1] = U_t
                    h_enc_sup[n+1] = U_t
                end
            else
                f8 = Interval.(factorial.(0:n_derivatives))
                t_left = interval(inf_t, 1)
                t_right = interval(1, sup_t)
                if n>=1
                    h_enc_inf[n+1] = 𝐅₈((t_left+1)/2, U_t_left, δ, Φ, para; δ₋)[n] / @interval(2)^n;
                else
                    h_enc_inf[n+1] = U_t_left
                end
                h_enc_sup[n+1] = (sum((derivatives_right[n+1:end] .* ((t_right - 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
            end
        elseif inf_t > 1
            f8 = Interval.(factorial.(0:8))
            h_enc_inf[n+1] = (sum((derivatives_right[n+1:end] .* ((t - 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
            h_enc_sup[n+1] = (sum((derivatives_right[n+1:end] .* ((t - 1) .^ (0:n_derivatives-n))) .* (f8[1:end-n] .^ (-1 * ones(n_derivatives+1 - n)))))
        end
    end
    h_enc = union.(h_enc_inf,h_enc_sup)
    return h_enc
end

function evaluation_hybrid_enclosure(t,U,δ,Φ::Sequence{CartesianPower{Chebyshev}},para,endpoint_derivatives;n_derivatives=8,δ₋=1)
    if inf(t)<=-1
        U_t_left = interval(zero(eltype(U(0))))
        if sup(t)<=-1
            U_t_right = interval(zero(eltype(U(0))))
        else
            t_right = interval(-1,sup(t))
            U_t_right = U(t_right)[1]
        end
    end
    if sup(t)>=1
        U_t_right = interval(zero(eltype(U(0))))
        if inf(t)>=1
            U_t_left = interval(zero(eltype(U(0))))
        else
            t_left = interval(inf(t),1)
            U_t_left = U(t_left)[1]
        end
    end
    if inf(t)>=-1 && sup(t)<=1
        U_t_left = interval(zero(eltype(U(0))))
        U_t_right = interval(zero(eltype(U(0))))
    end
    U_t = U(t)[1]
    h = hybrid_enclosure(t,U_t,U_t_left,U_t_right,δ,convert_Chebyshev_to_callable_function(Φ),para,endpoint_derivatives;n_derivatives,δ₋)
    return h
end

function evaluation_hybrid_enclosure(t,U,δ,Φ,para,endpoint_derivatives;n_derivatives=8,δ₋=1)
    if inf(t)<=-1
        U_t_left = interval(zero(eltype(U(0))))
        if sup(t)<=-1
            U_t_right = interval(zero(eltype(U(0))))
        else
            t_right = interval(-1,sup(t))
            U_t_right = U(t_right)[1]
        end
    end
    if sup(t)>=1
        U_t_right = interval(zero(eltype(U(0))))
        if inf(t)>=1
            U_t_left = interval(zero(eltype(U(0))))
        else
            t_left = interval(inf(t),1)
            U_t_left = U(t_left)[1]
        end
    end
    if inf(t)>=-1 && sup(t)<=1
        U_t_left = interval(zero(eltype(U(0))))
        U_t_right = interval(zero(eltype(U(0))))
    end
    U_t = U(t)[1]
    h = hybrid_enclosure(t,U_t,U_t_left,U_t_right,δ,Φ,para,endpoint_derivatives;n_derivatives,δ₋)
    return h
end

function Lipschitz_constant_f₁(t,u,δ,η₁,η₂,all_Φ::Sequence{CartesianPower{Chebyshev}},para)
    Iη₁ = interval(-1,1)*η₁;
    Iη₂ = interval(-1,1)*η₂;
    γ,κ,α,c = para;
    Φ,dΦ = eachcomponent(component(all_Φ,1:2));
    term1 = abs(δ + Iη₁)*abs(-γ + κ*c*dΦ(t*(δ+Iη₂) - α - c*(u+Iη₁))[0] );
    term2 = abs(γ*(u+Iη₁) + κ*Φ(t*(δ+Iη₂) - α - c*(u+Iη₁))[0]) + abs((δ + Iη₁)*κ*t*dΦ(t*(δ+Iη₂) - α - c*(u+Iη₁))[0]);
    return abs(Iη₁*term1 + Iη₂*term2)
end

function Lipschitz_constant_f₁(t,u,δ,η₁,η₂,all_Φ,para)
    Iη₁ = interval(-1,1)*η₁;
    Iη₂ = interval(-1,1)*η₂;
    γ,κ,α,c = para;
    all_Φ_eval = all_Φ(t*(δ+Iη₂) - α - c*(u+Iη₁));
    Φ,dΦ = all_Φ_eval[1:2];
    term1 = abs(δ + Iη₂)*abs(-γ + κ*c*dΦ );
    term2 = abs(γ*(u+Iη₁) + κ*Φ) + abs((δ + Iη₂)*κ*t*dΦ);
    return abs(Iη₁*term1 + Iη₂*term2)
end

function interpolate_solution_tight_D012(U_jℓ,δ,s_jℓ,Φ,C0_error,δ_error,N_new::Integer,Ni_new::T where T<:Real,
    𝐤,para,ppi::T where T<:Real;δ₋=1,sdiv=10,max_N=60,check_large_coeffs=1)
    # Interpolate ̄u and the bootstrapped derivative, directly. 
    # Outputs:
    # U_center : numerical center of u and its derivatives, as a Chebyshev()² object. 
    # U_radius : rigorous C⁰ error bounds for u and its derivative.
    # U : numerical center plus interval error bounds.
    𝐫₀ = interval(C0_error);  𝐫₁_sup = zero(eltype(U_jℓ));
    m = num_components(U_jℓ)
    im_inv = 1/@interval(m);    isdiv_inv = 1/@interval(sdiv)
    V = compute_V(U_jℓ,s_jℓ,Φ,δ,para,𝐤;δ₋=δ₋,sdiv=sdiv)
    Λ = 1 + (2 \ ppi) * log(Ni_new + 1)
    bound_scale = ppi*(1:𝐤).*(Ni_new .- (1:𝐤)).^(1:𝐤);
    while minimum(sup.(V * 4 ./ bound_scale))>sup(Λ*𝐫₀) && N_new<max_N
        N_new = N_new+1;
        Ni_new = Ni_new + 1;
        Λ = 1 + (2 \ ppi) * log(Ni_new + 1)
        bound_scale = ppi*(1:𝐤).*(Ni_new .- (1:𝐤)).^(1:𝐤);
    end
    k_best = findall(x->x==minimum(sup.(V * 4 ./ bound_scale)), sup.(V * 4 ./ bound_scale))[1]
    U, _, _, N_new, Ni_new = piecewise_to_smooth(U_jℓ,s_jℓ,N_new,Ni_new,ppi)
    if check_large_coeffs==1
        while sup(abs(coefficients(U)[end]))>1E-5 && N_new<max_N  # Increment N if there is a "large" high-order mode.
            U, _, _, N_new, Ni_new = piecewise_to_smooth(U_jℓ,s_jℓ,N_new+1,Ni_new+1,ppi;U_error = 0)
        end
    end
    _, _, U = interpolate_derivatives(U,δ,Φ,para,N_new,Ni_new,1,ppi;δ₋=δ₋)
    η₁ = C0_error;    η₂ = δ_error
    for j=1:m
        t0 = (j-1)*im_inv
        for k=1:sdiv
            t = t0 + interval(k-1,k)*isdiv_inv*im_inv
            t_shift = 2*interval(k-1,k)*isdiv_inv - 1
            u = component(U_jℓ,j)(t_shift)[0]
            𝐫₁_sup = max(𝐫₁_sup, sup( Lipschitz_constant_f₁(t,u,δ,η₁,η₂,Φ,para) ))
        end
    end
    U_center = copy(U)
    U_radius = interval.(zeros(eltype(U_jℓ),2))
    U_radius[1] = V[k_best] * 4 / (ppi * k_best * (Ni_new - k_best)^k_best) + Λ*𝐫₀
    U_radius[2] = V[k_best] * 4 / (ppi * k_best * (Ni_new - k_best + 1)^(k_best-1)) + Λ*interval(𝐫₁_sup)
    component(U,1)[0] += interval(-1,1)*U_radius[1]
    component(U,2)[0] += interval(-1,1)*U_radius[2]
    println("   ⋅ Interpolation completed successfully with N = $N_new.")
    return U_center, U_radius, U
end