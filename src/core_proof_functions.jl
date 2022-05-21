## Evaluation of G.

function 𝔉(u_jℓ,U_jℓ,δ,Φ,para, ppi::T where T<:Real; δ₋ = 1, check_edges=0, rigorous_error_control=0) # Interpolation of 𝐅ₚ.
    bstrap = 1;
    deg_Φ = length(Φ)-1;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = m*one(eltype(U_jℓ));
    degrees = compute_max_degree(k+1,deg_Φ);
    deg = degrees[bstrap];   degi = one(eltype(U_jℓ))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate 𝐅 on the domains ⋃ⱼ[s_j,s_{j+1}]
    components_𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    Fₚ_eval = zeros(eltype(U_jℓ),deg+3);
    if rigorous_error_control==1
        Fₚ_eval_error = Interval.(zeros(eltype(U_jℓ),m));
    else
        Fₚ_eval_error = zeros(eltype(U_jℓ),m);
    end
    for j ∈ 1:m
        for ℓ ∈ 1:deg+3
            Fₚ_eval[ℓ] = 𝐅( (j-1)/mi + (1/mi)*(nodes[ℓ]+1)/2, component(U_jℓ,j)(nodes[ℓ])[0] ,δ ,Φ ,para ; δ₋, check_edges)[1] 
        end
        if rigorous_error_control==1
            components_𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(Fₚ_eval)),deg+2,degi+2,ppi)
            Fₚ_eval_error[j] = Interval(-1,1)*norm(radius.(Fₚ_eval),Inf)
        else
            components_𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Fₚ_eval,deg+2,degi+2,ppi)
        end
    end
    # Rigorous error.
    Λ = 1 + (2/ppi)*log(deg+1); defect = Λ*Fₚ_eval_error;
    𝐅ₚ = Sequence(Chebyshev(deg+2)^m,components_𝐅ₚ) + Sequence(Chebyshev(0)^m,defect);
    return 𝐅ₚ
end

function 𝔉₈(u_jℓ,U_jℓ,δ,Φ,para, ppi::T where T<:Real, ind; δ₋=1, components_Φ=8) # Interpolation of 𝐅₈.
    deg_Φ = length(component(Φ,1))-1;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = m*one(eltype(U_jℓ));
    degrees = compute_max_degree(k+1,deg_Φ);
    deg = degrees[ind];   degi = one(eltype(U_jℓ))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate 𝐅 on the domains ⋃ⱼ[s_j,s_{j+1}]
    components_𝐅ₚ = Interval.(zeros(eltype(U_jℓ),(deg+3)*m));
    Fₚ_eval = Interval.(zeros(eltype(U_jℓ),deg+3));
    u = similar(component(U_jℓ,1));
    for j ∈ 1:m
        u[:] = coefficients(component(U_jℓ,j));
        for ℓ ∈ 1:deg+3
            Fₚ_eval[ℓ] = 𝐅₈( (j/mi)*(nodes[ℓ]+1)/2, u(nodes[ℓ])[0] ,δ ,Φ ,para,ind; δ₋, components_Φ)[1]  # 
        end
        components_𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Fₚ_eval,deg+2,degi+2,ppi)
    end
    𝐅ₚ = Sequence(Chebyshev(deg)^m,components_𝐅ₚ);
    return 𝐅ₚ
end

function 𝐅_limitLeft(u_jℓ,U_jℓ,δ,Φ,para; δ₋=1)
    # Evaluation of 𝐅 at all points sⱼ, with left limit for u.
    bstrap = 1;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = one(eltype(U_jℓ))*m;
    grid_F = zeros(eltype(U_jℓ),bstrap+1,m);
    for j ∈ 2:m
        grid_F[2:bstrap+1,j] .= 𝐅( (j-1)/mi, u_jℓ[k+2,j-1],δ ,Φ ,para; δ₋); # Note: the time argument is (j-1)/m.
        grid_F[1,j] = u_jℓ[k+2,j-1];
    end
    grid_F[2:bstrap+1,1] .= 𝐅( 0 , Φ(1)[0] ,δ ,Φ ,para; δ₋);
    grid_F[1,1] = Φ(1)[0];
    return grid_F
end

@views function gₚ(u_jℓ,U_jℓ,s_jℓ,δ,Φ,para, ppi::T where T<:Real;δ₋=1, check_edges=0,rigorous_error_control=0)
    bstrap = 1;
    k = length(component(U_jℓ,1))-2;
    ki = one(eltype(U_jℓ))*k;
    m = size(u_jℓ,2);
    mi = one(eltype(U_jℓ))*m;   Δ = 1/mi;
    fac_p = one(eltype(U_jℓ))*factorial.(0:bstrap-1);
    # Build the summation (pointwise eval)
    F_left = 𝐅_limitLeft(u_jℓ,U_jℓ,δ,Φ,para; δ₋)
    tsum = zeros(eltype(U_jℓ),k+2,m)
    for j∈1:m
        for q∈0:bstrap-1
            tsum[:,j] = tsum[:,j] + ((1/fac_p[q+1])*(s_jℓ[:,j].-s_jℓ[1,j]).^q)*F_left[q+1,j];
        end
    end
    # Build the integral (pointwise eval)
    Fₚ = 𝔉(u_jℓ,U_jℓ,δ,Φ,para,ppi; δ₋, check_edges,rigorous_error_control)
    tint = zeros(eltype(U_jℓ),k+2,m)
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    ∫ = RadiiPolynomial.Integral(1);
    igrand = similar(component(Fₚ,1));
    for j∈1:m
        igrand[:] = coefficients(component(Fₚ,j));
        intgrl = ∫*igrand;
        for ℓ∈2:k+2
            tint[ℓ,j] = ((Δ/2)^bstrap)/fac_p[bstrap]*(intgrl(nodes[ℓ]) - intgrl(-1))[0];
        end
    end
    # Return sum of summation and integral terms
    u_jℓ_new = tsum + tint;
    # Interpolate the result
    U_jℓ_new = convert_matrix_to_interp(u_jℓ_new,ppi);
    return u_jℓ_new , U_jℓ_new
end

function Gₚ(u_jℓ,U_jℓ,s_jℓ,δ,Φ,para, ppi::T where T<:Real; δ₋=1, check_edges=0,rigorous_error_control=0)
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);
    (g_jℓ,g_jℓ_interp) = gₚ(u_jℓ,U_jℓ,s_jℓ,δ,Φ,para,ppi; δ₋, check_edges,rigorous_error_control)
    h = δ - τ(0,u_jℓ[end,end],para);
    Gₚ_tuple = (g_jℓ - u_jℓ ,g_jℓ_interp - U_jℓ,h);
    Gₚ_vect = zeros(eltype(U_jℓ),m*(k+2)+1)
    Gₚ_vect[1:end-1] = reshape(g_jℓ - u_jℓ,m*(k+2));
    Gₚ_vect[end] = h;
    return Gₚ_tuple, Gₚ_vect
end

## Evaluation of derivatives.

function 𝐃𝐅(t,u,δ,Φ,para,ind; δ₋=1,components_Φ=2)
    # Evaluation of derivatives of bootstrapped vector fields, delta-scale, for specified argument t∈[0,1], u(t), given initial data Φ and parameter data para.
    Φ_arg = t*δ .- (para[3] + para[4]*u);         # argument of Φ in [-δ₋,0]
    Φ_arg = 1 .+ 2/δ₋*Φ_arg;                # scale argument to [-1,1]
    scale = ((2/δ₋) .^ (0:components_Φ-1))*one(eltype(Φ));    # compute derivative scale factor due to domain shift [-δ₋,0] to [-1,1]  
    𝐃𝐅 = DF(t,u,δ,scale.*Φ(Φ_arg),para); # Compute and δ-scale
end

function 𝔇𝔉(u_jℓ,U_jℓ,δ,Φ,para, ppi::T where T<:Real; δ₋=1,components_Φ=2,rigorous_error_control=0)  
    bstrap = 1;
    deg_Φ = length(component(Φ,1))-1;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = m*one(eltype(U_jℓ));
    degrees = compute_max_degree(k+1,deg_Φ);
    deg = compute_max_degree_DF(degrees[bstrap],k+1,deg_Φ);
    degi = one(eltype(U_jℓ))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate 𝐃𝐅 on the domains ⋃ⱼ[s_j,s_{j+1}]
    if rigorous_error_control==1
        D₁Fₚ_eval_error = Interval.(zeros(eltype(U_jℓ),m));
        D₂Fₚ_eval_error = Interval.(zeros(eltype(U_jℓ),m));
    else
        D₁Fₚ_eval_error = zeros(eltype(U_jℓ),m);
        D₂Fₚ_eval_error = zeros(eltype(U_jℓ),m);
    end
    components_D₁𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    components_D₂𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    D₁Fₚ_eval = zeros(eltype(U_jℓ),deg+3); D₂Fₚ_eval = zeros(eltype(U_jℓ),deg+3);
    for j ∈ 1:m
        for ℓ₁ ∈ 1:deg+3 
            @inbounds D₁Fₚ_eval[ℓ₁], D₂Fₚ_eval[ℓ₁] = 𝐃𝐅( (j-1)/mi + (1/mi)*(nodes[ℓ₁]+1)/2, component(U_jℓ,j)(nodes[ℓ₁])[0] ,δ ,Φ ,para,1; δ₋, components_Φ) 
        end
        if rigorous_error_control==1
            @inbounds components_D₁𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(D₁Fₚ_eval)),deg+2,degi+2,ppi)
            @inbounds components_D₂𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(D₂Fₚ_eval)),deg+2,degi+2,ppi)
            @inbounds D₁Fₚ_eval_error[j] = Interval(-1,1)*norm(radius.(D₁Fₚ_eval),Inf)
            @inbounds D₂Fₚ_eval_error[j] = Interval(-1,1)*norm(radius.(D₂Fₚ_eval),Inf)
        else
            @inbounds components_D₁𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(D₁Fₚ_eval,deg+2,degi+2,ppi)
            @inbounds components_D₂𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(D₂Fₚ_eval,deg+2,degi+2,ppi)
        end
    end
    # Rigorous error.
    Λ = 1 + (2/ppi)*log(deg+1); defect₁ = Λ*D₁Fₚ_eval_error;    defect₂ = Λ*D₂Fₚ_eval_error;

    𝐃₁𝐅ₚ = Sequence(Chebyshev(deg+2)^m,components_D₁𝐅ₚ) + Sequence(Chebyshev(0)^m,defect₁);
    𝐃₂𝐅ₚ = Sequence(Chebyshev(deg+2)^m,components_D₂𝐅ₚ) + Sequence(Chebyshev(0)^m,defect₂);
    if bstrap==1
        return [𝐃₁𝐅ₚ], [𝐃₂𝐅ₚ]
    end
end

function DG_sum_and_BC(DF_interp,u_jℓ,U_jℓ,δ,s_jℓ,m,k,Φ,para)
    bstrap = 1;
    𝐃₁𝐅, 𝐃₂𝐅 = DF_interp;
    Dsum = zeros(eltype(s_jℓ),m*(k+2)+1,m*(k+2)+1)
    fac_p = one(eltype(s_jℓ))*factorial.(0:bstrap-1);
    # Derivatives with respect to u(s_jℓ); all rows except boundarycondition.
    for j ∈ 1:m-1
        for ℓ∈1:k+2
            Dsum[1+(j)*(k+2) + ℓ-1,(k+2)*j] = one(eltype(s_jℓ));
            for q∈1:bstrap-1
                Dsum[1+(j)*(k+2) + ℓ-1,(k+2)*j] += ((1/fac_p[q+1])*(s_jℓ[ℓ,j]-s_jℓ[1,j])^q)*component(𝐃₁𝐅[q],j)(1)[0];
            end
        end
    end
    # Derivative with respect to δ; all rows except boundarycondition.
    for j∈1:m
        for ℓ∈1:k+2
            for q∈1:bstrap-1
                Dsum[1+(j-1)*(k+2)+(ℓ-1),1+m*(k+2)] += ((1/fac_p[q+1])*(s_jℓ[ℓ,j]-s_jℓ[1,j])^q)*component(𝐃₂𝐅[q],j)(1)[0];
            end
        end
    end
    # boundarycondition
    Dsum[end,end-1] = one(eltype(s_jℓ))*(-D₂τ(1,component(U_jℓ,m)(1)[0],para));
    Dsum[end,end] = one(eltype(s_jℓ));
    return Dsum
end

function DG_integral(DF_interp,u_jℓ,U_jℓ,δ,s_jℓ,m,k,Φ,para,ppi::T where T<:Real)
    bstrap = 1;
    𝐃₁𝐅, 𝐃₂𝐅 = DF_interp;   one! = one(eltype(U_jℓ));
    Dint = zeros(eltype(s_jℓ),m*(k+2)+1,m*(k+2)+1)*one!;
    fac_p =  one(eltype(s_jℓ))*factorial(bstrap-1);
    s_id = Sequence(Chebyshev(1),one(eltype(s_jℓ))*[0;1/2]);
    ki = one(eltype(s_jℓ))*k;   mi = one(eltype(s_jℓ))*m;   Δ = 1/mi;
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    ∫ = RadiiPolynomial.Integral(1);
    # Derivatives w.r.t. u(s_jℓ).
    if bstrap==1
        igrand = similar(component(𝐃₁𝐅[1],1));
    elseif bstrap==2
        igrand = similar(component(𝐃₁𝐅[1],1)*s_id^(bstrap-1));
    end
    e_vect = zeros(eltype(s_jℓ),k+2);
    Cheb_matrix = Cheb_node_to_spectral_matrix(e_vect,nodes)
    for ℓ∈1:k+2
        for j∈1:m
            if bstrap==1
                igrand[:] = coefficients(component(𝐃₁𝐅[1],j));
            elseif bstrap==2
                pow_s = ((nodes[ℓ] - s_id)^(bstrap-1))
                igrand[:] = coefficients(pow_s*component(𝐃₁𝐅[bstrap],j));
            end
            for i∈1:k+2
                e_vect[:] = I[1:k+2,i];
                e_cheb = component(convert_matrix_to_interp(e_vect,Cheb_matrix),1);
                intgrl = ∫*(igrand*e_cheb);
                Dint[(j-1)*(k+2) + ℓ,(j-1)*(k+2)+i] = ((Δ/2)^bstrap)/fac_p*(intgrl(nodes[ℓ]) - intgrl(-1))[0];
            end
        end
    end 
    # Derivatives w.r.t δ
    if bstrap==1
        igrand = similar(component(𝐃₂𝐅[1],1));
    elseif bstrap==2
        igrand = similar(component(𝐃₂𝐅[1],1)*s_id^(bstrap-1));
    end
    for ℓ∈1:k+2
        for j∈1:m
            if bstrap==1
                igrand[:] = coefficients(component(𝐃₂𝐅[1],j));
            elseif bstrap==2
                pow_s = ((nodes[ℓ] - s_id)^(bstrap-1))
                igrand[:] = coefficients(pow_s*component(𝐃₂𝐅[bstrap],j));
            end
            intgrl = ∫*(igrand);
            Dint[(j-1)*(k+2) + ℓ, end] = ((Δ/2)^bstrap)/fac_p*(intgrl(nodes[ℓ]) - intgrl(-1))[0];
        end
    end
    return Dint
end

function DGₚ(u_tuple,δ,s_jℓ,Φ,para, ppi::T where T<:Real; δ₋=1,components_Φ=2,rigorous_error_control=0)
    u_jℓ = u_tuple[1];  U_jℓ = u_tuple[2];
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);
    D = 𝔇𝔉(u_jℓ,U_jℓ,δ,Φ,para,ppi;δ₋,components_Φ,rigorous_error_control);
    Dsum = DG_sum_and_BC(D, u_jℓ, U_jℓ, δ, s_jℓ, m, k, Φ, para);
    Dint = DG_integral(D, u_jℓ, U_jℓ, δ, s_jℓ, m, k, Φ, para, ppi);
    Id = one(eltype(s_jℓ))*I[1:m*(k+2)+1,1:m*(k+2)+1];    Id[end,end]=zero(eltype(s_jℓ));
    return Dsum + Dint - Id
end

## Newton.

function Newton(u_tuple,δ,s_jℓ,Φ,para, ppi::T where T<:Real,tol,maxiter; δ₋=1)
    u_jℓ,U_jℓ = u_tuple;
    x = [reshape(u_jℓ,size(u_jℓ,1)*size(u_jℓ,2),1); δ];
    G,G_vect = Gₚ(u_jℓ,U_jℓ,s_jℓ,δ,component(Φ,1),para,ppi; δ₋);
    DG = DGₚ((Float64.(u_jℓ),Float64.(U_jℓ)),Float64.(δ),Float64.(s_jℓ),Float64.(component(Φ,1:2)),Float64.(para),Float64.(ppi); δ₋);
    iter = 1;
    correction = DG\G_vect;
    print("\r   ⋅ Newton iteration 1: residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(DG*G_vect)),sigdigits=2) ); println();
    while norm(correction,1)>tol && iter < maxiter
        print("\r   ⋅ Newton iteration ", iter+1, ": "); 
        x[:] = x - correction;
        u_jℓ[:] = reshape(x[1:end-1],size(u_jℓ,1),size(u_jℓ,2));
        U_jℓ = convert_matrix_to_interp(u_jℓ,ppi);
        δ = x[end];
        G,G_vect = Gₚ(u_jℓ,U_jℓ,s_jℓ,δ,component(Φ,1),para,ppi;δ₋);
        DG = DGₚ((Float64.(u_jℓ),Float64.(U_jℓ)),Float64.(δ),Float64.(s_jℓ),Float64.(component(Φ,1:2)),Float64.(para),Float64.(ppi);δ₋);
        correction = DG\G_vect;
        print("residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(correction,1)),sigdigits=2) ); println()
        iter += 1;
    end
    residual = norm(G_vect)
    return u_jℓ,U_jℓ,δ,residual
end

function Newton_multiple_Phi(u_tuple,δ,s_jℓ,Φ,Φ_low,para, ppi::T where T<:Real,tol,maxiter; δ₋=1)
    # Newton's method running with multiple representations of Φ, for efficiency in subsequent proofs.
    # Computes F using Φ, computes DF using Φ_low.
    u_jℓ,U_jℓ = u_tuple;
    x = [reshape(u_jℓ,size(u_jℓ,1)*size(u_jℓ,2),1); δ];
    G,G_vect = Gₚ(u_jℓ,U_jℓ,s_jℓ,δ,component(Φ,1),para,ppi; δ₋);
    DG = DGₚ((Float64.(u_jℓ),Float64.(U_jℓ)),Float64.(δ),Float64.(s_jℓ),Float64.(component(Φ_low,1:2)),Float64.(para),Float64.(ppi); δ₋);
    iter = 1;
    correction = DG\G_vect;
    print("\r   ⋅ Newton iteration 1: residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(DG*G_vect)),sigdigits=2) ); println();
    while norm(correction,1)>tol && iter < maxiter
        print("\r   ⋅ Newton iteration ", iter+1, ": "); 
        x[:] = x - correction;
        u_jℓ[:] = reshape(x[1:end-1],size(u_jℓ,1),size(u_jℓ,2));
        U_jℓ = convert_matrix_to_interp(u_jℓ,ppi);
        δ = x[end];
        G,G_vect = Gₚ(u_jℓ,U_jℓ,s_jℓ,δ,component(Φ,1),para,ppi;δ₋);
        DG = DGₚ((Float64.(u_jℓ),Float64.(U_jℓ)),Float64.(δ),Float64.(s_jℓ),Float64.(component(Φ_low,1:2)),Float64.(para),Float64.(ppi); δ₋);
        correction = DG\G_vect;
        print("residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(correction,1)),sigdigits=2) ); println()
        iter += 1;
    end
    return u_jℓ,U_jℓ,δ
end

## ---- Computer-assisted proofs; bounds ----
# Norm.
function norm_ℓν¹(v)
    norm([one(eltype(v));2*ones(eltype(v),length(v)-1)].*v,1);
end

# Vector norm
function prenorm_vec(g,r₀)
    weight = [ones(eltype(g),length(g)-1); (1/r₀)];
    pnv = (abs.(g)).*weight;
    return pnv
end

# Operator norm
function prenorm_op(M,r₀)
    pno = zeros(eltype(M),size(M,1));
    weight = [ones(eltype(M),size(M,1)-1); r₀];
    for j=1:size(M,1)-1
        pno[j] = dot(abs.(M[j,:]),weight);
    end
    pno[end] = (1/r₀)*dot(abs.(M[end,:]),weight);
    return pno
end

# Lag interpolation
function interpolate_lag(u_tuple,δ,para)
    u_jℓ,U_jℓ = u_tuple;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = @interval m;
    if k==0
        error("Only k>0 is supported.")
    end
    lagⱼ = similar(component(U_jℓ,1));
    lag_interp_coeffs = zeros(eltype(U_jℓ),(k+2)*m);
    for j∈1:m
        tδ = Sequence(Chebyshev(1),[(j-1)*δ/mi;0] + (δ/2)*(1/mi)*[1;1/2]);
        lagⱼ[:] = coefficients(tδ - τ(component(U_jℓ,j),para));
        lag_interp_coeffs[1+(j-1)*(k+2):j*(k+2)] = coefficients(lagⱼ);
    end
    lag_interp = Sequence(Chebyshev(k+1)^m,lag_interp_coeffs);
    return lag_interp
end

# Second derivative interpolation
function 𝐃²𝐅(t,u,δ,Φ,para;δ₋=1)
    # Evaluation of second derivatives of bootstrapped vector fields, delta-scale, for specified argument t∈[0,1], u(t), given initial data Φ and parameter data para.
    Φ_arg = t*δ .- (para[3] + para[4]*u);   # argument of Φ in [-δ₋,0]
    Φ_arg = 1 .+ (2/δ₋)*Φ_arg;                   # scale argument to [-1,1]
    scale = ((2/δ₋) .^ (0:3))*one(eltype(Φ));    # compute derivative scale factor due to domain shift [-δ₋,0] to [-1,1]   
    𝐃²𝐅 = D²F(t,u,δ,scale.*Φ(Φ_arg),para); # Compute and δ-scale
end

function 𝔇²𝔉(u_jℓ,U_jℓ,δ,Φ,para, ppi::T where T<:Real; δ₋=1)   # Interpolation of D²𝐅ₚ.
    bstrap == 1;
    deg_Φ = length(component(Φ,1))-1;
    k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = m*one(eltype(U_jℓ));
    degrees = compute_max_degree(k+1,deg_Φ);
    deg = compute_max_degree_D²F(degrees[bstrap],k+1,deg_Φ);
    degi = one(eltype(U_jℓ))*deg;
    nodes = cheb_nodes(deg,degi,ppi);
    # Interpolate 𝐃𝐅 on the domains ⋃ⱼ[s_j,s_{j+1}]
    components_Dx²𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    components_Dxδ𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    components_Dδ²𝐅ₚ = zeros(eltype(U_jℓ),(deg+3)*m);
    Dx²Fₚ_eval = zeros(eltype(U_jℓ),deg+3); DxδFₚ_eval = zeros(eltype(U_jℓ),deg+3); 
    Dδ²Fₚ_eval = zeros(eltype(U_jℓ),deg+3);
    u = similar(component(U_jℓ,1));
    for j ∈ 1:m
        u[:] = coefficients(component(U_jℓ,j));
        for ℓ₁ ∈ 1:deg+3
            Dx²Fₚ_eval[ℓ₁], DxδFₚ_eval[ℓ₁], Dδ²Fₚ_eval[ℓ₁] = 𝐃²𝐅( (j-1)/mi + (1/mi)*(nodes[ℓ₁]+1)/2, component(U_jℓ,j)(nodes[ℓ₁])[0] ,δ ,Φ ,para; δ₋)
        end
        components_Dx²𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Dx²Fₚ_eval,deg,degi,ppi);
        components_Dxδ𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DxδFₚ_eval,deg,degi,ppi);
        components_Dδ²𝐅ₚ[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Dδ²Fₚ_eval,deg,degi,ppi);
    end
    𝐃x²𝐅ₚ = Sequence(Chebyshev(deg)^m,components_Dx²𝐅ₚ);
    𝐃xδ𝐅ₚ = Sequence(Chebyshev(deg)^m,components_Dxδ𝐅ₚ);
    𝐃δ²𝐅ₚ = Sequence(Chebyshev(deg)^m,components_Dδ²𝐅ₚ);
    return [𝐃x²𝐅ₚ 𝐃xδ𝐅ₚ 𝐃δ²𝐅ₚ] 
end

# Check lag inclusion of numerical solution
function lag_range(u_tuple,δ,para,error;δ₋=1)
    lag = interpolate_lag(u_tuple,δ,para);
    u_jℓ,U_jℓ = u_tuple;    m = size(u_jℓ,2);
    ercheb = similar(lag(zero(eltype(U_jℓ))));  ercheb[:] .= para[4]*error;
    lag_eval = lag(Interval(-1,1)) - ercheb;
    range_lag = Interval(-1,1)*ones(Interval,m);
    for j∈1:m
        maxv = sup(lag_eval[j]); minv = inf(lag_eval[j]);
        max_translate = (2/δ₋)*Interval(maxv)+1;
        min_translate = (2/δ₋)*Interval(minv)+1;
        range_lag[j] = Interval(inf(min_translate),sup(max_translate));
    end
    return range_lag, lag
end

# Approximate inverse
function compute_A(DG,T::Type)
    A = interval.(inv(T.(mid.(DG))));
    return A
end

# Y bound
function bound_Y(A,G_vect,r₀;do_round=0)
    if do_round==0
        Y = prenorm_vec(A*G_vect,r₀);
    elseif do_round==1
        G! = interval.(Float64.(inf.(G_vect),RoundDown),Float64.(sup.(G_vect),RoundUp));
        Y = prenorm_vec(A*G!,r₀);
    end
    return Y
end

# Y∞ bound
@views function bound_Y∞(u_tuple,δ,Φ::Sequence{CartesianPower{Chebyshev}},para;δ₋=1,sdiv=10)   
    bstrap = 1;
    u_jℓ, U_jℓ = u_tuple;   k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = Interval(m);   Δ = 1/mi;
    bound_dᵏF = zeros(eltype(U_jℓ),m);
    U = Array{Sequence}(undef,k+1);  U[1] = U_jℓ;
    U_eval = zeros(eltype(U_jℓ),k+1)
    scale_dΦ = ((2/δ₋) .^ (0:k))*one(eltype(Φ));     #derivative scaling for Φ
    scale_dU = (2*m) .^ (0:k)*one(eltype(U_jℓ));     #derivative scaling for U
    for n=1:k+1-bstrap
        U[n+1] = Derivative(1)*U[n];
    end
    isd = interval(sdiv)
    for j∈1:m
        for jj=1:sdiv
            tⱼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tⱼ_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Φ_arg = δ*tⱼ - para[3] - para[4]*component(U[1],j)(tⱼ_cheb)[0];
            Φ_arg = 1 .+ 2/δ₋*Φ_arg;    #rescale to correct domain [-1,1]
            for e=1:k+1
                U_eval[e] = component(U[e],j)(tⱼ_cheb)[0];
            end
            bound_dᵏF[j] = interval( max( sup(bound_dᵏF[j]), sup(abs(F₈(scale_dU.*U_eval[1:k+1],scale_dΦ.*coefficients(component(Φ,1:k+1)(Φ_arg)),para,δ;order=k)[k+1])) ) );
        end
    end
    Cₙ = 1/(Interval(factorial(k+1))*2^(2*(k+bstrap-1)));
    Y∞ = Δ^(k+bstrap)*Cₙ*norm(sup.(bound_dᵏF),Inf);
    return Interval(sup(Y∞))
end

@views function bound_Y∞(u_tuple,δ,Φ,para;δ₋=1,sdiv=10)   
    bstrap = 1;
    u_jℓ, U_jℓ = u_tuple;   k = length(component(U_jℓ,1))-2;
    m = size(u_jℓ,2);   mi = Interval(m);   Δ = 1/mi;
    bound_dᵏF = zeros(eltype(U_jℓ),m);
    U = Array{Sequence}(undef,k+1);  U[1] = U_jℓ;
    U_eval = zeros(eltype(U_jℓ),k+1)
    scale_dΦ = ((2/δ₋) .^ (0:k))*one(eltype( Φ(interval(0))[1] ));     #derivative scaling for Φ
    scale_dU = (2*m) .^ (0:k)*one(eltype(U_jℓ));     #derivative scaling for U
    for n=1:k+1-bstrap
        U[n+1] = Derivative(1)*U[n];
    end
    isd = interval(sdiv)
    for j∈1:m
        for jj=1:sdiv
            tⱼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tⱼ_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Φ_arg = δ*tⱼ - para[3] - para[4]*component(U[1],j)(tⱼ_cheb)[0];
            Φ_arg = 1 .+ 2/δ₋*Φ_arg;    #rescale to correct domain [-1,1]
            for e=1:k+1
                U_eval[e] = component(U[e],j)(tⱼ_cheb)[0];
            end
            bound_dᵏF[j] = interval( max( sup(bound_dᵏF[j]), sup(abs(F₈(scale_dU.*U_eval[1:k+1],scale_dΦ.*Φ(Φ_arg)[1:k+1],para,δ;order=k)[k+1])) ) );
        end
    end
    Cₙ = 1/(Interval(factorial(k+1))*2^(2*(k+bstrap-1)));
    Y∞ = Δ^(k+bstrap)*Cₙ*norm(sup.(bound_dᵏF),Inf);
    return Interval(sup(Y∞))
end

# Z0 bound
function bound_Z0(A,DG,r₀;do_round=0)
    n = size(A)[1]
    wrap_A = LinearOperator(ParameterSpace()^n,ParameterSpace()^n,A)    # Wrap A and DG (later) as LinearOperators for fast multiplication.
    if do_round==0
        wrap_DG = LinearOperator(ParameterSpace()^n,ParameterSpace()^n,DG)
        defect = I - wrap_A*wrap_DG;
    elseif do_round==1
        DG! = interval.(Float64.(inf.(DG),RoundDown),Float64.(sup.(DG),RoundUp));
        wrap_DG! = LinearOperator(ParameterSpace()^n,ParameterSpace()^n,DG!)
        defect = I - wrap_A*wrap_DG!
    end
    Z0 = prenorm_op(defect.coefficients,r₀);
    return Z0
end

# Z1 bound
function bound_Z1(A,u_tuple,δ,s_jℓ,Φ,para,ppi,r₀,r_∞; δ₋=1,sdiv=10)
    bstrap = 1;
    u_jℓ, U_jℓ = u_tuple;   k = size(u_jℓ,1)-2;    ki = Interval(k);
    m = size(u_jℓ,2);
    ρ = zeros(eltype(u_jℓ),k+2,m);
    fac_p = one(eltype(u_jℓ))*factorial(bstrap);
    for ℓ∈1:k+2
        for j∈1:m
            t = range(inf(s_jℓ[ℓ,j]) , sup(s_jℓ[1,j]) ;length=sdiv);
            DF_bound = zero(eltype(u_jℓ));
            for sd = 1:sdiv-1
                t_sd = @interval(t[sd])∪@interval(t[sd+1]);
                DF_bound = max( DF_bound, abs( 𝐃𝐅(t_sd,component(U_jℓ,j)(-1+2*t_sd)[0], δ, Φ,para,1 ;δ₋ )[1] ) )
            end
            ρ[ℓ,j] = (abs(s_jℓ[ℓ,1] - 0)^bstrap)/fac_p * DF_bound;
        end
    end
    ρ_vect = reshape(ρ,(k+2)*(m),1);
    Aρ = [abs.(view(A,1:(k+2)*m,1:(k+2)*m))*ρ_vect;0];
    return r_∞*prenorm_vec(Aρ,r₀)
end

# Z2 bound
function bound_Z2(A,u_tuple,δ,s_jℓ,Φ::Sequence{CartesianPower{Chebyshev}},para,ppi,r,r₀,r_∞;δ₋=1)
    bstrap = 1;
    u_jℓ, _ = u_tuple;   k = size(u_jℓ,1)-2;    m = size(u_jℓ,2); ki = Interval(k);
    Λ = 1 + (2\ppi)*log(ki+1);  mi = interval(m);
    ϱ₂ = Interval(0);   # Note: the hessian of the lag function vanishes FOR OUR EXAMPLE. This line of code is NOT GENERAL.
    ϱ₁ = zeros(eltype(u_jℓ),k+2,m); u_rad = r*(Λ+r_∞);
    Φ_bound = zeros(eltype(u_jℓ),3);    δ_rad = r*r₀*Interval(-1,1);
    range_lag,_ = lag_range(u_tuple,δ+δ_rad,para,u_rad*Interval(-1,1); δ₋);   ## MOVE THIS COMPUTATION OUTSIDE OF THE Z2 BOUND.
    # Get inclusion for Φ
    for n∈1:3
        for j∈1:m
            Φ_bound[n] = Interval(max(sup(Φ_bound[n]),sup(component(Φ,n)(range_lag[j])[0])));
        end
    end
    bound_2derivative = zeros(eltype(u_jℓ),3,bstrap);
    for q∈1:bstrap
        for j=1:m
            u_enc = Interval(-1,1)*(Λ*norm(u_jℓ[:,j],Inf) + u_rad);
            bound_2derivative[:,q] .= max.( bound_2derivative[:,q], D²F(Interval(inf(s_jℓ[1,j]),sup(s_jℓ[end,j])),u_enc,δ+δ_rad,Interval(-1,1)*Φ_bound,para) );
        end
    end
    Hess_Θ = zeros(eltype(u_jℓ),2,2);
    for q∈0:bstrap-1
        Hess_Θ[:] = abs.([bound_2derivative[1,q+1] bound_2derivative[2,q+1]; 
                    bound_2derivative[2,q+1] bound_2derivative[3,q+1] ]);
    end
    Δs_jℓ = similar(s_jℓ);
    for ℓ∈1:k+2
        Δs_jℓ[ℓ,:] = s_jℓ[ℓ,:] - s_jℓ[1,:]
    end
    int_multi = dot(Interval(-1,1)*[Λ+r_∞;r₀],Hess_Θ*(Interval(-1,1)*[Λ+r_∞;r₀]))/factorial(bstrap)
    int_bound = Δs_jℓ.^bstrap*int_multi;
    ϱ₁ = reshape(int_bound,m*(k+2));
    ϱ = [ϱ₁;ϱ₂];
    Z2 = prenorm_vec(A*ϱ*r,r₀);
    return Z2
end

function bound_Z2(A,u_tuple,δ,s_jℓ,Φ,para,ppi,r,r₀,r_∞;δ₋=1)
    bstrap = 1;
    u_jℓ, _ = u_tuple;   k = size(u_jℓ,1)-2;    m = size(u_jℓ,2); ki = Interval(k);
    Λ = 1 + (2\ppi)*log(ki+1);  mi = interval(m);
    ϱ₂ = Interval(0);   # Note: the hessian of the lag function vanishes FOR OUR EXAMPLE. This line of code is NOT GENERAL.
    ϱ₁ = zeros(eltype(u_jℓ),k+2,m); u_rad = r*(Λ+r_∞);
    Φ_bound = zeros(eltype(u_jℓ),3);    δ_rad = r*r₀*Interval(-1,1);
    range_lag,_ = lag_range(u_tuple,δ+δ_rad,para,u_rad*Interval(-1,1); δ₋);   ## MOVE THIS COMPUTATION OUTSIDE OF THE Z2 BOUND.
    # Get inclusion for Φ
    for n∈1:3
        for j∈1:m
            Φ_bound[n] = Interval(max(sup(Φ_bound[n]),Φ(range_lag[j])[n]));
        end
    end
    bound_2derivative = zeros(eltype(u_jℓ),3,bstrap);
    for q∈1:bstrap
        for j=1:m
            u_enc = Interval(-1,1)*(Λ*norm(u_jℓ[:,j],Inf) + u_rad);
            bound_2derivative[:,q] .= max.( bound_2derivative[:,q], D²F(Interval(inf(s_jℓ[1,j]),sup(s_jℓ[end,j])),u_enc,δ+δ_rad,Interval(-1,1)*Φ_bound,para) );
        end
    end
    Hess_Θ = zeros(eltype(u_jℓ),2,2);
    for q∈0:bstrap-1
        Hess_Θ[:] = abs.([bound_2derivative[1,q+1] bound_2derivative[2,q+1]; 
                    bound_2derivative[2,q+1] bound_2derivative[3,q+1] ]);
    end
    Δs_jℓ = similar(s_jℓ);
    for ℓ∈1:k+2
        Δs_jℓ[ℓ,:] = s_jℓ[ℓ,:] - s_jℓ[1,:]
    end
    int_multi = dot(Interval(-1,1)*[Λ+r_∞;r₀],Hess_Θ*(Interval(-1,1)*[Λ+r_∞;r₀]))/factorial(bstrap)
    int_bound = Δs_jℓ.^bstrap*int_multi;
    ϱ₁ = reshape(int_bound,m*(k+2));
    ϱ = [ϱ₁;ϱ₂];
    Z2 = prenorm_vec(A*ϱ*r,r₀);
    return Z2
end

# Z∞ bound
function bound_Z∞(u_tuple,δ,Φ::Sequence{CartesianPower{Chebyshev}},para,ppi,r,r₀,r_∞; δ₋=1, sdiv=10)
    bstrap = 1;
    u_jℓ, U_jℓ = u_tuple;   k = size(u_jℓ,1)-2;    ki = Interval(k);
    Λ = 1 + (2\ppi)*log(ki+1);
    Cₖₚ_star_1 = (1+Λ)*(ppi/4)^bstrap*factorial(k+1-bstrap)/Interval(factorial(k+1));
    Cₖₚ_star_2 = 1/(Interval(factorial(bstrap))^2^bstrap)   # Note: always have k≥p and p≤2.
    Cₖₚ = Interval(min(sup(Cₖₚ_star_1),sup(Cₖₚ_star_2)));
    m = size(u_jℓ,2);    mi = Interval(m);
    δ_rad = r*r₀*Interval(-1,1);    u_rad = r*(Λ+r_∞);
    u_perturb = u_rad*Interval(-1,1);
    bounds_DΘ = zeros(eltype(U_jℓ),m);
    scale_dΦ = ((2/δ₋) .^ (0:1))*one(eltype(Φ));  
    isd = interval(sdiv);
    for j∈1:m
        for jj=1:sdiv
            tⱼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tⱼ_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Φ_arg = δ*tⱼ - para[3] - para[4]*component(U_jℓ,j)(tⱼ_cheb)[0];
            Φ_arg = 1 .+ 2/δ₋*Φ_arg;    #rescale to correct domain [-1,1]
            U_eval = component(U_jℓ,j)(tⱼ_cheb)[0];
            𝔇𝔉1, 𝔇𝔉2 = DF(tⱼ,U_eval + u_perturb,δ + δ_rad,scale_dΦ.*coefficients(Φ(Φ_arg)),para);
            bounds_DΘ[j] = abs(𝔇𝔉1)*(Λ+r_∞) + abs(𝔇𝔉2)*(Λ+r_∞);
        end
    end
    Z∞ = Cₖₚ*(1/mi)^bstrap*norm(bounds_DΘ,Inf);
    return Z∞
end

function bound_Z∞(u_tuple,δ,Φ,para,ppi,r,r₀,r_∞; δ₋=1, sdiv=10)
    bstrap = 1;
    u_jℓ, U_jℓ = u_tuple;   k = size(u_jℓ,1)-2;    ki = Interval(k);
    Λ = 1 + (2\ppi)*log(ki+1);
    Cₖₚ_star_1 = (1+Λ)*(ppi/4)^bstrap*factorial(k+1-bstrap)/Interval(factorial(k+1));
    Cₖₚ_star_2 = 1/(Interval(factorial(bstrap))^2^bstrap) 
    Cₖₚ = Interval(min(sup(Cₖₚ_star_1),sup(Cₖₚ_star_2)));
    m = size(u_jℓ,2);    mi = Interval(m);
    δ_rad = r*r₀*Interval(-1,1);    u_rad = r*(Λ+r_∞);
    u_perturb = u_rad*Interval(-1,1);
    bounds_DΘ = zeros(eltype(U_jℓ),m);
    scale_dΦ = ((2/δ₋) .^ (0:1))*one(eltype(Φ(interval(0))[1]));  
    isd = interval(sdiv);
    for j∈1:m
        for jj=1:sdiv
            tⱼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tⱼ_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Φ_arg = δ*tⱼ - para[3] - para[4]*component(U_jℓ,j)(tⱼ_cheb)[0];
            Φ_arg = 1 .+ 2/δ₋*Φ_arg;    #rescale to correct domain [-1,1]
            U_eval = component(U_jℓ,j)(tⱼ_cheb)[0];
            𝔇𝔉1, 𝔇𝔉2 = DF(tⱼ,U_eval + u_perturb,δ + δ_rad,scale_dΦ.*Φ(Φ_arg)[1:2],para);
            bounds_DΘ[j] = abs(𝔇𝔉1)*(Λ+r_∞) + abs(𝔇𝔉2)*(Λ+r_∞);
        end
    end
    Z∞ = Cₖₚ*(1/mi)^bstrap*norm(bounds_DΘ,Inf);
    return Z∞
end

# radii polynomial
function radpol(G,A,DG,u,δ,Φ,para,s_jℓ,ppi,r_star,r₀,r_∞;δ₋=1,rounding_YZ=0,debug=0,sdiv=10)
    if rounding_YZ==1
        u_jℓ! = interval.(Float64.(inf.(u[1]),RoundDown),Float64.(sup.(u[1]),RoundUp));
        U_jℓ! = Sequence(space(u[2]),interval.(Float64.(inf.(coefficients(u[2])),RoundDown),Float64.(sup.(coefficients(u[2])),RoundUp)));
        u! = (u_jℓ!,U_jℓ!);
        δ! = interval(Float64(inf(δ),RoundDown),Float64(sup(δ),RoundUp))
        Φ! = Sequence(space(Φ),interval.(Float64.(inf.(coefficients(Φ)),RoundDown),Float64.(sup.(coefficients(Φ)),RoundUp)));
        s_jℓ! = interval.(Float64.(inf.(s_jℓ),RoundDown),Float64.(sup.(s_jℓ),RoundUp));
        ppi! = interval(Float64(inf(ppi),RoundDown),Float64(sup(ppi),RoundUp))
        para! = interval.(Float64.(inf.(para),RoundDown),Float64.(sup.(para),RoundUp));
        δ₋! = interval(Float64(inf(δ₋),RoundDown),Float64(sup(δ₋),RoundUp))
        r_star! = interval(Float64(inf(r_star),RoundDown),Float64(sup(r_star),RoundUp))
        r₀! = interval(Float64(inf(r₀),RoundDown),Float64(sup(r₀),RoundUp))
        r_∞! = interval(Float64(inf(r_∞),RoundDown),Float64(sup(r_∞),RoundUp))
    elseif rounding_YZ==0
        u! = u; δ! = δ; Φ! = Φ; s_jℓ! = s_jℓ;   ppi! = ppi; para! = para;   δ₋! = δ₋;   r_star! = r_star;   r₀! = r₀;   r_∞! = r_∞;
    end
    if debug==0
        Y0 = bound_Y(A,G,r₀!;do_round=rounding_YZ);   sY = norm(sup.(Y0),Inf)
        println("   ⋅ Y₀ bound computed: $sY.")
        Y∞ = bound_Y∞(u!,δ!,Φ!,para!;δ₋,sdiv);   sY∞ = sup(Y∞);
        println("   ⋅ Y∞ bound computed: $sY∞.")
        Z0 = bound_Z0(A,DG,r₀!;do_round=rounding_YZ); sZ0 = norm(sup.(Z0),Inf);
        println("   ⋅ Z₀ bound computed: $sZ0.")
        Z1 = bound_Z1(A,u!,δ!,s_jℓ!,component(Φ!,1:2),para!,ppi!,r₀!,r_∞!;δ₋=δ₋!,sdiv);  sZ1 = norm(sup.(Z1),Inf);
        println("   ⋅ Z₁ bound computed: $sZ1.")
        Z2 = bound_Z2(A,u!,δ!,s_jℓ!,component(Φ!,1:3),para!,ppi!,r_star!,r₀!,r_∞!;δ₋=δ₋!);   sZ2 = norm(sup.(Z2),Inf);
        println("   ⋅ Z₂ bound computed: $sZ2.")
        Z∞ = bound_Z∞(u!,δ!,component(Φ!,1:2),para!,ppi!,r_star!,r₀!,r_∞!;δ₋=δ₋!,sdiv);  sZ∞ = sup(Z∞);
        println("   ⋅ Z∞ bound computed: $sZ∞.")
        poly_finite = [Y0 (Z2+Z1+Z0 .- 1)];
        poly_∞ = [norm(Y∞,Inf) Z∞ - r_∞];
        r = [- poly_finite[:,1]./poly_finite[:,2] ; -poly_∞[1]/poly_∞[2] ];
        if minimum(inf.(r))<0
            println("Existence interval is empty; a radius is negative.")
            ie = ∅
            C0_err = Inf;   δ_err = Inf;
            # @infiltrate
        elseif maximum(nextfloat.(sup.(r)))>r_star
            println("Existence interval is empty; a radius is larger than r_star.")
            ie = ∅
            C0_err = Inf;   δ_err = Inf;
            # @infiltrate
        else
            ie = Interval(maximum(nextfloat.(sup.(r))),inf(r_star))
            k = size(u[1],1)-2;    ki = Interval(k);
            Λ = 1 + (2\ppi)*log(ki+1);
            C0_err = inf(ie)*(Λ+r_∞); δ_err = inf(ie)*r₀;   sC0 = sup(C0_err);  sδ = sup(δ_err);
            println("   ⋅ C⁰ enclosure: $sC0.")
            println("   ⋅ δ enlosure: $sδ.")
        end
        return Y0, Y∞, Z0, Z1, Z2, Z∞, poly_finite, poly_∞, r, ie, C0_err, δ_err
    else
        # @infiltrate
        return []
    end
end

function radpol(G,A,DG,u,δ,Φ_Cheb,Φ_function,para,s_jℓ,ppi,r_star,r₀,r_∞;δ₋=1,rounding_YZ=0,debug=0,sdiv=10)
    if rounding_YZ==1
        u_jℓ! = interval.(Float64.(inf.(u[1]),RoundDown),Float64.(sup.(u[1]),RoundUp));
        U_jℓ! = Sequence(space(u[2]),interval.(Float64.(inf.(coefficients(u[2])),RoundDown),Float64.(sup.(coefficients(u[2])),RoundUp)));
        u! = (u_jℓ!,U_jℓ!);
        δ! = interval(Float64(inf(δ),RoundDown),Float64(sup(δ),RoundUp))
        Φ_Cheb! = Sequence(space(Φ_Cheb),interval.(Float64.(inf.(coefficients(Φ_Cheb)),RoundDown),Float64.(sup.(coefficients(Φ_Cheb)),RoundUp)));
        s_jℓ! = interval.(Float64.(inf.(s_jℓ),RoundDown),Float64.(sup.(s_jℓ),RoundUp));
        ppi! = interval(Float64(inf(ppi),RoundDown),Float64(sup(ppi),RoundUp))
        para! = interval.(Float64.(inf.(para),RoundDown),Float64.(sup.(para),RoundUp));
        δ₋! = interval(Float64(inf(δ₋),RoundDown),Float64(sup(δ₋),RoundUp))
        r_star! = interval(Float64(inf(r_star),RoundDown),Float64(sup(r_star),RoundUp))
        r₀! = interval(Float64(inf(r₀),RoundDown),Float64(sup(r₀),RoundUp))
        r_∞! = interval(Float64(inf(r_∞),RoundDown),Float64(sup(r_∞),RoundUp))
        Φ_function! = t-> interval.(Float64.(inf.(Φ_function(t)),RoundDown),Float64.(sup.(Φ_function(t)),RoundUp));
    elseif rounding_YZ==0
        u! = u; δ! = δ; Φ_Cheb! = Φ_Cheb; s_jℓ! = s_jℓ;   ppi! = ppi; para! = para;   δ₋! = δ₋; r_star! = r_star;   r₀! = r₀;   r_∞! = r_∞; Φ_function! = Φ_function;
    end
    if debug==0
        Y0 = bound_Y(A,G,r₀!;do_round=rounding_YZ);   sY = norm(sup.(Y0),Inf)
        println("   ⋅ Y₀ bound computed: $sY.")
        Y∞ = bound_Y∞(u!,δ!,Φ_function!,para!;δ₋,sdiv);   sY∞ = sup(Y∞);
        println("   ⋅ Y∞ bound computed: $sY∞.")
        Z0 = bound_Z0(A,DG,r₀;do_round=rounding_YZ); sZ0 = norm(sup.(Z0),Inf);
        println("   ⋅ Z₀ bound computed: $sZ0.")
        Z1 = bound_Z1(A,u!,δ!,s_jℓ!,component(Φ_Cheb!,1:2),para!,ppi!,r₀!,r_∞!;δ₋=δ₋!,sdiv);  sZ1 = norm(sup.(Z1),Inf);
        println("   ⋅ Z₁ bound computed: $sZ1.")
        Z2 = bound_Z2(A,u!,δ!,s_jℓ!,Φ_function!,para!,ppi!,r_star!,r₀!,r_∞!;δ₋=δ₋!);   sZ2 = norm(sup.(Z2),Inf);
        println("   ⋅ Z₂ bound computed: $sZ2.")
        Z∞ = bound_Z∞(u!,δ!,Φ_function!,para!,ppi!,r_star!,r₀!,r_∞!;δ₋=δ₋!,sdiv);  sZ∞ = sup(Z∞);
        println("   ⋅ Z∞ bound computed: $sZ∞.")
        poly_finite = [Y0 (Z2+Z1+Z0 .- 1)];
        poly_∞ = [norm(Y∞,Inf) Z∞ - r_∞];
        r = [- poly_finite[:,1]./poly_finite[:,2] ; -poly_∞[1]/poly_∞[2] ];
        if minimum(inf.(r))<0
            println("Existence interval is empty; a radius is negative.")
            ie = ∅
            C0_err = Inf;   δ_err = Inf;
            # @infiltrate
        elseif maximum(nextfloat.(sup.(r)))>r_star
            println("Existence interval is empty; a radius is larger than r_star.")
            ie = ∅
            C0_err = Inf;   δ_err = Inf;
            # infiltrate
        else
            ie = Interval(maximum(nextfloat.(sup.(r))),inf(r_star))
            k = size(u[1],1)-2;    ki = Interval(k);
            Λ = 1 + (2\ppi)*log(ki+1);
            C0_err = inf(ie)*(Λ+r_∞); δ_err = inf(ie)*r₀;   sC0 = sup(C0_err);  sδ = sup(δ_err);
            println("   ⋅ C⁰ enclosure: $sC0.")
            println("   ⋅ δ enlosure: $sδ.")
        end
        return Y0, Y∞, Z0, Z1, Z2, Z∞, poly_finite, poly_∞, r, ie, C0_err, δ_err
    else
        # @infiltrate
        return []
    end
end

function modify_candidate_zero(iu,δ,para,δ₋;desired_tolerance=1e-20)
    α = para[3];    c = para[4]
    u_0 = iu[1][1,1];    u_1 = iu[1][end,end]
    scale = sup.( [δ -c*u_1 ; 0 -c*u_0]\[-desired_tolerance - (δ - α - c*u_1) ; -δ₋ + desired_tolerance - (-α - c*u_0)] )
    while sup((1+scale[1])*δ - α - (1+scale[2])*c*u_1)>0 && inf((1/δ₋)*(-α - (1+scale[2])*c*u_0))<-1
        scale = scale*1.01;
    end
    iu_mod = ( interval(1+scale[2])*iu[1], Sequence(space(iu[2]),interval(1+scale[2])*coefficients(iu[2])) );
    δ_mod = (1+scale[2])*δ;
    return iu_mod,δ_mod
end

function modify_candidate_zero_rightside_only(iu,δ,para;desired_tolerance=1e-20)
    α = para[3];    c = para[4]
    u_1 = iu[1][end,end]
    if sup(δ - α - c*u_1)>0
        ϵ₁ = (-interval(desired_tolerance) - (δ - α - c*u_1))/δ;
        while sup((1+ϵ₁)*δ - α -c*u_1)>0
            ϵ₁ = ϵ₁*1.01;
        end
        δ_mod = (1+ϵ₁)*δ
    else
        δ_mod = δ
    end
    return iu,δ_mod
end

function check_lag_monotonicity(Φ::Sequence{CartesianPower{Chebyshev}},U,s_jℓ,δ,r_C0)
    # Note 1, this code is example-specific.
    γ,κ,α,c = para
    t = LinRange(-1,1,sdiv)
    d_lag = Inf
    du = interval.(zero(eltype(U)))
    u = interval.(zero(eltype(U)))
    Φ_eval = interval.(zero(eltype(U)))
    dΦ_eval = interval.(zero(eltype(U)))
    for n=1:num_components(U)
        u = component(U,n)(interval(-1,1))[0] + r_C0*interval(-1,1)
        t = interval(inf(s_jℓ[1,n]),sup(s_jℓ[end,n]))
        Φ_eval = component(Φ,1)(t*δ - α - c*u)[0]
        dΦ_eval = component(Φ,2)(t*δ - α - c*u)[0]
        du = δ*(-γ*u - κ*Φ_eval)
        d_lag = min(d_lag, inf(δ - c*du) )
    end
    check = d_lag>0
    println("   ⋅ The lag is monotone: $check.")
    return d_lag>0
end

function check_lag_monotonicity(para,Φ,U,s_jℓ,δ,r_C0)
    # Note 1, this code is example-specific.
    γ,κ,α,c = para
    d_lag = Inf
    du = interval.(zero(eltype(U)))
    u = interval.(zero(eltype(U)))
    Φ_eval = interval.(zero(eltype(U)))
    dΦ_eval = interval.(zero(eltype(U)))
    for n=1:num_components(U)
        u = component(U,n)(interval(-1,1))[0] + r_C0*interval(-1,1)
        t = interval(inf(s_jℓ[1,n]),sup(s_jℓ[end,n]))
        Φ_eval = Φ(t*δ - α - c*u)[1]
        dΦ_eval = Φ(t*δ - α - c*u)[2]
        du = δ*(-γ*u - κ*Φ_eval)
        d_lag = min(d_lag, inf(δ - c*du) )
    end
    check = d_lag>0
    println("   ⋅ The lag is monotone: $check.")
    return d_lag>0
end

function proof_section_radpol_interp_monotonicity(G_n,A_n,DG_n,iu_n,iδ_n,Φ,Φ_function,ipara,s_jℓ_n,ppi,r_star,r_0,r_∞,δ_val_previous,kΦ,kΦ_low,𝐤,rounding_YZ)
    println(": Starting evaluation of the bounds, radii polynomials. :")
    _, _, _, _, _, _, _, _, _, _, C0_err, δ_err = radpol(G_n,A_n,DG_n,iu_n,iδ_n,Φ,Φ_function,ipara,s_jℓ_n,ppi,r_star,r_0,r_∞;δ₋=δ_val_previous,rounding_YZ);
    println(": Starting high-accuracy interpolation. :")
    _,_,Φ_interp_high_accuracy = interpolate_solution_tight_D012(iu_n[2],iδ_n,s_jℓ_n,Φ_function,C0_err,δ_err,kΦ,interval(kΦ),𝐤,ipara,ppi;δ₋=δ_val_previous,sdiv=10,check_large_coeffs=1)
    println(": Starting low-accuracy interpolation. :")
    _,_,Φ_interp_low_accuracy = interpolate_solution_tight_D012(iu_n[2],iδ_n,s_jℓ_n,Φ_function,C0_err,δ_err,kΦ_low,interval(kΦ_low),𝐤,ipara,ppi;δ₋=δ_val_previous,sdiv=10,max_N=kΦ_low,check_large_coeffs=0)
    println(": Building hybrid enclosure. :")
    Φ_function_next = t -> evaluation_hybrid_enclosure(t,convert_Chebyshev_to_callable_function(component(Φ_interp_high_accuracy,1)),
                            iδ_n + interval(-1,1)*δ_err,Φ_function,ipara,all_derivatives(iu_n[1],ipara,Φ_function,iδ_n + interval(-1,1)*δ_err;
                            δ₋=δ_val_previous);n_derivatives=𝐤,δ₋=δ_val_previous)
    println(": Verifying monotonicity of time lag. :")
    check_lag_monotonicity(ipara,Φ_function,iu_n[2],s_jℓ_n,iδ_n + interval(-1,1)*δ_err,C0_err)
    return C0_err, δ_err, Φ_interp_high_accuracy, Φ_interp_low_accuracy, Φ_function_next
end

function proof_section_radpol_interp_monotonicity_nohybrid(G_n,A_n,DG_n,iu_n,iδ_n,Φ,Φ_function,ipara,s_jℓ_n,ppi,r_star,r_0,r_∞,δ_val_previous,kΦ,kΦ_low,𝐤,rounding_YZ)
    println(": Starting evaluation of the bounds, radii polynomials. :")
    _, _, _, _, _, _, _, _, _, _, C0_err, δ_err = radpol(G_n,A_n,DG_n,iu_n,iδ_n,Φ,Φ_function,ipara,s_jℓ_n,ppi,r_star,r_0,r_∞;δ₋=δ_val_previous,rounding_YZ);
    println(": Starting high-accuracy interpolation. :")
    _,_,Φ_interp_high_accuracy = interpolate_solution_tight_D012(iu_n[2],iδ_n,s_jℓ_n,Φ_function,C0_err,δ_err,kΦ,interval(kΦ),𝐤,ipara,ppi;δ₋=δ_val_previous,sdiv=10,check_large_coeffs=1)
    println(": Starting low-accuracy interpolation. :")
    _,_,Φ_interp_low_accuracy = interpolate_solution_tight_D012(iu_n[2],iδ_n,s_jℓ_n,Φ_function,C0_err,δ_err,kΦ_low,interval(kΦ_low),𝐤,ipara,ppi;δ₋=δ_val_previous,sdiv=10,max_N=kΦ_low,check_large_coeffs=0)
    println(": Verifying monotonicity of time lag. :")
    check_lag_monotonicity(ipara,Φ_function,iu_n[2],s_jℓ_n,iδ_n + interval(-1,1)*δ_err,C0_err)
    return C0_err, δ_err, Φ_interp_high_accuracy, Φ_interp_low_accuracy
end
