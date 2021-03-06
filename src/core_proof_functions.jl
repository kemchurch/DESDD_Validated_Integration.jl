## Evaluation of G.

function š(u_jā,U_jā,Ī“,Ī¦,para, ppi::T where T<:Real; Ī“ā = 1, check_edges=0, rigorous_error_control=0) # Interpolation of šā.
    bstrap = 1;
    deg_Ī¦ = length(Ī¦)-1;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = m*one(eltype(U_jā));
    degrees = compute_max_degree(k+1,deg_Ī¦);
    deg = degrees[bstrap];   degi = one(eltype(U_jā))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate š on the domains āā±¼[s_j,s_{j+1}]
    components_šā = zeros(eltype(U_jā),(deg+3)*m);
    Fā_eval = zeros(eltype(U_jā),deg+3);
    if rigorous_error_control==1
        Fā_eval_error = Interval.(zeros(eltype(U_jā),m));
    else
        Fā_eval_error = zeros(eltype(U_jā),m);
    end
    for j ā 1:m
        for ā ā 1:deg+3
            Fā_eval[ā] = š( (j-1)/mi + (1/mi)*(nodes[ā]+1)/2, component(U_jā,j)(nodes[ā])[0] ,Ī“ ,Ī¦ ,para ; Ī“ā, check_edges)[1] 
        end
        if rigorous_error_control==1
            components_šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(Fā_eval)),deg+2,degi+2,ppi)
            Fā_eval_error[j] = Interval(-1,1)*norm(radius.(Fā_eval),Inf)
        else
            components_šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Fā_eval,deg+2,degi+2,ppi)
        end
    end
    # Rigorous error.
    Ī = 1 + (2/ppi)*log(deg+1); defect = Ī*Fā_eval_error;
    šā = Sequence(Chebyshev(deg+2)^m,components_šā) + Sequence(Chebyshev(0)^m,defect);
    return šā
end

function šā(u_jā,U_jā,Ī“,Ī¦,para, ppi::T where T<:Real, ind; Ī“ā=1, components_Ī¦=8) # Interpolation of šā.
    deg_Ī¦ = length(component(Ī¦,1))-1;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = m*one(eltype(U_jā));
    degrees = compute_max_degree(k+1,deg_Ī¦);
    deg = degrees[ind];   degi = one(eltype(U_jā))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate š on the domains āā±¼[s_j,s_{j+1}]
    components_šā = Interval.(zeros(eltype(U_jā),(deg+3)*m));
    Fā_eval = Interval.(zeros(eltype(U_jā),deg+3));
    u = similar(component(U_jā,1));
    for j ā 1:m
        u[:] = coefficients(component(U_jā,j));
        for ā ā 1:deg+3
            Fā_eval[ā] = šā( (j/mi)*(nodes[ā]+1)/2, u(nodes[ā])[0] ,Ī“ ,Ī¦ ,para,ind; Ī“ā, components_Ī¦)[1]  # 
        end
        components_šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Fā_eval,deg+2,degi+2,ppi)
    end
    šā = Sequence(Chebyshev(deg)^m,components_šā);
    return šā
end

function š_limitLeft(u_jā,U_jā,Ī“,Ī¦,para; Ī“ā=1)
    # Evaluation of š at all points sā±¼, with left limit for u.
    bstrap = 1;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = one(eltype(U_jā))*m;
    grid_F = zeros(eltype(U_jā),bstrap+1,m);
    for j ā 2:m
        grid_F[2:bstrap+1,j] .= š( (j-1)/mi, u_jā[k+2,j-1],Ī“ ,Ī¦ ,para; Ī“ā); # Note: the time argument is (j-1)/m.
        grid_F[1,j] = u_jā[k+2,j-1];
    end
    grid_F[2:bstrap+1,1] .= š( 0 , Ī¦(1)[0] ,Ī“ ,Ī¦ ,para; Ī“ā);
    grid_F[1,1] = Ī¦(1)[0];
    return grid_F
end

@views function gā(u_jā,U_jā,s_jā,Ī“,Ī¦,para, ppi::T where T<:Real;Ī“ā=1, check_edges=0,rigorous_error_control=0)
    bstrap = 1;
    k = length(component(U_jā,1))-2;
    ki = one(eltype(U_jā))*k;
    m = size(u_jā,2);
    mi = one(eltype(U_jā))*m;   Ī = 1/mi;
    fac_p = one(eltype(U_jā))*factorial.(0:bstrap-1);
    # Build the summation (pointwise eval)
    F_left = š_limitLeft(u_jā,U_jā,Ī“,Ī¦,para; Ī“ā)
    tsum = zeros(eltype(U_jā),k+2,m)
    for jā1:m
        for qā0:bstrap-1
            tsum[:,j] = tsum[:,j] + ((1/fac_p[q+1])*(s_jā[:,j].-s_jā[1,j]).^q)*F_left[q+1,j];
        end
    end
    # Build the integral (pointwise eval)
    Fā = š(u_jā,U_jā,Ī“,Ī¦,para,ppi; Ī“ā, check_edges,rigorous_error_control)
    tint = zeros(eltype(U_jā),k+2,m)
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    ā« = RadiiPolynomial.Integral(1);
    igrand = similar(component(Fā,1));
    for jā1:m
        igrand[:] = coefficients(component(Fā,j));
        intgrl = ā«*igrand;
        for āā2:k+2
            tint[ā,j] = ((Ī/2)^bstrap)/fac_p[bstrap]*(intgrl(nodes[ā]) - intgrl(-1))[0];
        end
    end
    # Return sum of summation and integral terms
    u_jā_new = tsum + tint;
    # Interpolate the result
    U_jā_new = convert_matrix_to_interp(u_jā_new,ppi);
    return u_jā_new , U_jā_new
end

function Gā(u_jā,U_jā,s_jā,Ī“,Ī¦,para, ppi::T where T<:Real; Ī“ā=1, check_edges=0,rigorous_error_control=0)
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);
    (g_jā,g_jā_interp) = gā(u_jā,U_jā,s_jā,Ī“,Ī¦,para,ppi; Ī“ā, check_edges,rigorous_error_control)
    h = Ī“ - Ļ(0,u_jā[end,end],para);
    Gā_tuple = (g_jā - u_jā ,g_jā_interp - U_jā,h);
    Gā_vect = zeros(eltype(U_jā),m*(k+2)+1)
    Gā_vect[1:end-1] = reshape(g_jā - u_jā,m*(k+2));
    Gā_vect[end] = h;
    return Gā_tuple, Gā_vect
end

## Evaluation of derivatives.

function šš(t,u,Ī“,Ī¦,para,ind; Ī“ā=1,components_Ī¦=2)
    # Evaluation of derivatives of bootstrapped vector fields, delta-scale, for specified argument tā[0,1], u(t), given initial data Ī¦ and parameter data para.
    Ī¦_arg = t*Ī“ .- (para[3] + para[4]*u);         # argument of Ī¦ in [-Ī“ā,0]
    Ī¦_arg = 1 .+ 2/Ī“ā*Ī¦_arg;                # scale argument to [-1,1]
    scale = ((2/Ī“ā) .^ (0:components_Ī¦-1))*one(eltype(Ī¦));    # compute derivative scale factor due to domain shift [-Ī“ā,0] to [-1,1]  
    šš = DF(t,u,Ī“,scale.*Ī¦(Ī¦_arg),para); # Compute and Ī“-scale
end

function šš(u_jā,U_jā,Ī“,Ī¦,para, ppi::T where T<:Real; Ī“ā=1,components_Ī¦=2,rigorous_error_control=0)  
    bstrap = 1;
    deg_Ī¦ = length(component(Ī¦,1))-1;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = m*one(eltype(U_jā));
    degrees = compute_max_degree(k+1,deg_Ī¦);
    deg = compute_max_degree_DF(degrees[bstrap],k+1,deg_Ī¦);
    degi = one(eltype(U_jā))*deg;
    nodes = cheb_nodes(deg+2,degi+2,ppi);
    # Interpolate šš on the domains āā±¼[s_j,s_{j+1}]
    if rigorous_error_control==1
        DāFā_eval_error = Interval.(zeros(eltype(U_jā),m));
        DāFā_eval_error = Interval.(zeros(eltype(U_jā),m));
    else
        DāFā_eval_error = zeros(eltype(U_jā),m);
        DāFā_eval_error = zeros(eltype(U_jā),m);
    end
    components_Dāšā = zeros(eltype(U_jā),(deg+3)*m);
    components_Dāšā = zeros(eltype(U_jā),(deg+3)*m);
    DāFā_eval = zeros(eltype(U_jā),deg+3); DāFā_eval = zeros(eltype(U_jā),deg+3);
    for j ā 1:m
        for āā ā 1:deg+3 
            @inbounds DāFā_eval[āā], DāFā_eval[āā] = šš( (j-1)/mi + (1/mi)*(nodes[āā]+1)/2, component(U_jā,j)(nodes[āā])[0] ,Ī“ ,Ī¦ ,para,1; Ī“ā, components_Ī¦) 
        end
        if rigorous_error_control==1
            @inbounds components_Dāšā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(DāFā_eval)),deg+2,degi+2,ppi)
            @inbounds components_Dāšā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(Interval.(mid.(DāFā_eval)),deg+2,degi+2,ppi)
            @inbounds DāFā_eval_error[j] = Interval(-1,1)*norm(radius.(DāFā_eval),Inf)
            @inbounds DāFā_eval_error[j] = Interval(-1,1)*norm(radius.(DāFā_eval),Inf)
        else
            @inbounds components_Dāšā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DāFā_eval,deg+2,degi+2,ppi)
            @inbounds components_Dāšā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DāFā_eval,deg+2,degi+2,ppi)
        end
    end
    # Rigorous error.
    Ī = 1 + (2/ppi)*log(deg+1); defectā = Ī*DāFā_eval_error;    defectā = Ī*DāFā_eval_error;

    šāšā = Sequence(Chebyshev(deg+2)^m,components_Dāšā) + Sequence(Chebyshev(0)^m,defectā);
    šāšā = Sequence(Chebyshev(deg+2)^m,components_Dāšā) + Sequence(Chebyshev(0)^m,defectā);
    if bstrap==1
        return [šāšā], [šāšā]
    end
end

function DG_sum_and_BC(DF_interp,u_jā,U_jā,Ī“,s_jā,m,k,Ī¦,para)
    bstrap = 1;
    šāš, šāš = DF_interp;
    Dsum = zeros(eltype(s_jā),m*(k+2)+1,m*(k+2)+1)
    fac_p = one(eltype(s_jā))*factorial.(0:bstrap-1);
    # Derivatives with respect to u(s_jā); all rows except boundarycondition.
    for j ā 1:m-1
        for āā1:k+2
            Dsum[1+(j)*(k+2) + ā-1,(k+2)*j] = one(eltype(s_jā));
            for qā1:bstrap-1
                Dsum[1+(j)*(k+2) + ā-1,(k+2)*j] += ((1/fac_p[q+1])*(s_jā[ā,j]-s_jā[1,j])^q)*component(šāš[q],j)(1)[0];
            end
        end
    end
    # Derivative with respect to Ī“; all rows except boundarycondition.
    for jā1:m
        for āā1:k+2
            for qā1:bstrap-1
                Dsum[1+(j-1)*(k+2)+(ā-1),1+m*(k+2)] += ((1/fac_p[q+1])*(s_jā[ā,j]-s_jā[1,j])^q)*component(šāš[q],j)(1)[0];
            end
        end
    end
    # boundarycondition
    Dsum[end,end-1] = one(eltype(s_jā))*(-DāĻ(1,component(U_jā,m)(1)[0],para));
    Dsum[end,end] = one(eltype(s_jā));
    return Dsum
end

function DG_integral(DF_interp,u_jā,U_jā,Ī“,s_jā,m,k,Ī¦,para,ppi::T where T<:Real)
    bstrap = 1;
    šāš, šāš = DF_interp;   one! = one(eltype(U_jā));
    Dint = zeros(eltype(s_jā),m*(k+2)+1,m*(k+2)+1)*one!;
    fac_p =  one(eltype(s_jā))*factorial(bstrap-1);
    s_id = Sequence(Chebyshev(1),one(eltype(s_jā))*[0;1/2]);
    ki = one(eltype(s_jā))*k;   mi = one(eltype(s_jā))*m;   Ī = 1/mi;
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    ā« = RadiiPolynomial.Integral(1);
    # Derivatives w.r.t. u(s_jā).
    if bstrap==1
        igrand = similar(component(šāš[1],1));
    elseif bstrap==2
        igrand = similar(component(šāš[1],1)*s_id^(bstrap-1));
    end
    e_vect = zeros(eltype(s_jā),k+2);
    Cheb_matrix = Cheb_node_to_spectral_matrix(e_vect,nodes)
    for āā1:k+2
        for jā1:m
            if bstrap==1
                igrand[:] = coefficients(component(šāš[1],j));
            elseif bstrap==2
                pow_s = ((nodes[ā] - s_id)^(bstrap-1))
                igrand[:] = coefficients(pow_s*component(šāš[bstrap],j));
            end
            for iā1:k+2
                e_vect[:] = I[1:k+2,i];
                e_cheb = component(convert_matrix_to_interp(e_vect,Cheb_matrix),1);
                intgrl = ā«*(igrand*e_cheb);
                Dint[(j-1)*(k+2) + ā,(j-1)*(k+2)+i] = ((Ī/2)^bstrap)/fac_p*(intgrl(nodes[ā]) - intgrl(-1))[0];
            end
        end
    end 
    # Derivatives w.r.t Ī“
    if bstrap==1
        igrand = similar(component(šāš[1],1));
    elseif bstrap==2
        igrand = similar(component(šāš[1],1)*s_id^(bstrap-1));
    end
    for āā1:k+2
        for jā1:m
            if bstrap==1
                igrand[:] = coefficients(component(šāš[1],j));
            elseif bstrap==2
                pow_s = ((nodes[ā] - s_id)^(bstrap-1))
                igrand[:] = coefficients(pow_s*component(šāš[bstrap],j));
            end
            intgrl = ā«*(igrand);
            Dint[(j-1)*(k+2) + ā, end] = ((Ī/2)^bstrap)/fac_p*(intgrl(nodes[ā]) - intgrl(-1))[0];
        end
    end
    return Dint
end

function DGā(u_tuple,Ī“,s_jā,Ī¦,para, ppi::T where T<:Real; Ī“ā=1,components_Ī¦=2,rigorous_error_control=0)
    u_jā = u_tuple[1];  U_jā = u_tuple[2];
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);
    D = šš(u_jā,U_jā,Ī“,Ī¦,para,ppi;Ī“ā,components_Ī¦,rigorous_error_control);
    Dsum = DG_sum_and_BC(D, u_jā, U_jā, Ī“, s_jā, m, k, Ī¦, para);
    Dint = DG_integral(D, u_jā, U_jā, Ī“, s_jā, m, k, Ī¦, para, ppi);
    Id = one(eltype(s_jā))*I[1:m*(k+2)+1,1:m*(k+2)+1];    Id[end,end]=zero(eltype(s_jā));
    return Dsum + Dint - Id
end

## Newton.

function Newton(u_tuple,Ī“,s_jā,Ī¦,para, ppi::T where T<:Real,tol,maxiter; Ī“ā=1)
    u_jā,U_jā = u_tuple;
    x = [reshape(u_jā,size(u_jā,1)*size(u_jā,2),1); Ī“];
    G,G_vect = Gā(u_jā,U_jā,s_jā,Ī“,component(Ī¦,1),para,ppi; Ī“ā);
    DG = DGā((Float64.(u_jā),Float64.(U_jā)),Float64.(Ī“),Float64.(s_jā),Float64.(component(Ī¦,1:2)),Float64.(para),Float64.(ppi); Ī“ā);
    iter = 1;
    correction = DG\G_vect;
    print("\r   ā Newton iteration 1: residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(DG*G_vect)),sigdigits=2) ); println();
    while norm(correction,1)>tol && iter < maxiter
        print("\r   ā Newton iteration ", iter+1, ": "); 
        x[:] = x - correction;
        u_jā[:] = reshape(x[1:end-1],size(u_jā,1),size(u_jā,2));
        U_jā = convert_matrix_to_interp(u_jā,ppi);
        Ī“ = x[end];
        G,G_vect = Gā(u_jā,U_jā,s_jā,Ī“,component(Ī¦,1),para,ppi;Ī“ā);
        DG = DGā((Float64.(u_jā),Float64.(U_jā)),Float64.(Ī“),Float64.(s_jā),Float64.(component(Ī¦,1:2)),Float64.(para),Float64.(ppi);Ī“ā);
        correction = DG\G_vect;
        print("residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(correction,1)),sigdigits=2) ); println()
        iter += 1;
    end
    residual = norm(G_vect)
    return u_jā,U_jā,Ī“,residual
end

function Newton_multiple_Phi(u_tuple,Ī“,s_jā,Ī¦,Ī¦_low,para, ppi::T where T<:Real,tol,maxiter; Ī“ā=1)
    # Newton's method running with multiple representations of Ī¦, for efficiency in subsequent proofs.
    # Computes F using Ī¦, computes DF using Ī¦_low.
    u_jā,U_jā = u_tuple;
    x = [reshape(u_jā,size(u_jā,1)*size(u_jā,2),1); Ī“];
    G,G_vect = Gā(u_jā,U_jā,s_jā,Ī“,component(Ī¦,1),para,ppi; Ī“ā);
    DG = DGā((Float64.(u_jā),Float64.(U_jā)),Float64.(Ī“),Float64.(s_jā),Float64.(component(Ī¦_low,1:2)),Float64.(para),Float64.(ppi); Ī“ā);
    iter = 1;
    correction = DG\G_vect;
    print("\r   ā Newton iteration 1: residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(DG*G_vect)),sigdigits=2) ); println();
    while norm(correction,1)>tol && iter < maxiter
        print("\r   ā Newton iteration ", iter+1, ": "); 
        x[:] = x - correction;
        u_jā[:] = reshape(x[1:end-1],size(u_jā,1),size(u_jā,2));
        U_jā = convert_matrix_to_interp(u_jā,ppi);
        Ī“ = x[end];
        G,G_vect = Gā(u_jā,U_jā,s_jā,Ī“,component(Ī¦,1),para,ppi;Ī“ā);
        DG = DGā((Float64.(u_jā),Float64.(U_jā)),Float64.(Ī“),Float64.(s_jā),Float64.(component(Ī¦_low,1:2)),Float64.(para),Float64.(ppi); Ī“ā);
        correction = DG\G_vect;
        print("residual(F) = ", round(Float64(norm(G_vect)),sigdigits=2), ", residual(AF) = ", round(Float64(norm(correction,1)),sigdigits=2) ); println()
        iter += 1;
    end
    return u_jā,U_jā,Ī“
end

## ---- Computer-assisted proofs; bounds ----
# Norm.
function norm_āĪ½Ā¹(v)
    norm([one(eltype(v));2*ones(eltype(v),length(v)-1)].*v,1);
end

# Vector norm
function prenorm_vec(g,rā)
    weight = [ones(eltype(g),length(g)-1); (1/rā)];
    pnv = (abs.(g)).*weight;
    return pnv
end

# Operator norm
function prenorm_op(M,rā)
    pno = zeros(eltype(M),size(M,1));
    weight = [ones(eltype(M),size(M,1)-1); rā];
    for j=1:size(M,1)-1
        pno[j] = dot(abs.(M[j,:]),weight);
    end
    pno[end] = (1/rā)*dot(abs.(M[end,:]),weight);
    return pno
end

# Lag interpolation
function interpolate_lag(u_tuple,Ī“,para)
    u_jā,U_jā = u_tuple;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = @interval m;
    if k==0
        error("Only k>0 is supported.")
    end
    lagā±¼ = similar(component(U_jā,1));
    lag_interp_coeffs = zeros(eltype(U_jā),(k+2)*m);
    for jā1:m
        tĪ“ = Sequence(Chebyshev(1),[(j-1)*Ī“/mi;0] + (Ī“/2)*(1/mi)*[1;1/2]);
        lagā±¼[:] = coefficients(tĪ“ - Ļ(component(U_jā,j),para));
        lag_interp_coeffs[1+(j-1)*(k+2):j*(k+2)] = coefficients(lagā±¼);
    end
    lag_interp = Sequence(Chebyshev(k+1)^m,lag_interp_coeffs);
    return lag_interp
end

# Second derivative interpolation
function šĀ²š(t,u,Ī“,Ī¦,para;Ī“ā=1)
    # Evaluation of second derivatives of bootstrapped vector fields, delta-scale, for specified argument tā[0,1], u(t), given initial data Ī¦ and parameter data para.
    Ī¦_arg = t*Ī“ .- (para[3] + para[4]*u);   # argument of Ī¦ in [-Ī“ā,0]
    Ī¦_arg = 1 .+ (2/Ī“ā)*Ī¦_arg;                   # scale argument to [-1,1]
    scale = ((2/Ī“ā) .^ (0:3))*one(eltype(Ī¦));    # compute derivative scale factor due to domain shift [-Ī“ā,0] to [-1,1]   
    šĀ²š = DĀ²F(t,u,Ī“,scale.*Ī¦(Ī¦_arg),para); # Compute and Ī“-scale
end

function šĀ²š(u_jā,U_jā,Ī“,Ī¦,para, ppi::T where T<:Real; Ī“ā=1)   # Interpolation of DĀ²šā.
    bstrap == 1;
    deg_Ī¦ = length(component(Ī¦,1))-1;
    k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = m*one(eltype(U_jā));
    degrees = compute_max_degree(k+1,deg_Ī¦);
    deg = compute_max_degree_DĀ²F(degrees[bstrap],k+1,deg_Ī¦);
    degi = one(eltype(U_jā))*deg;
    nodes = cheb_nodes(deg,degi,ppi);
    # Interpolate šš on the domains āā±¼[s_j,s_{j+1}]
    components_DxĀ²šā = zeros(eltype(U_jā),(deg+3)*m);
    components_DxĪ“šā = zeros(eltype(U_jā),(deg+3)*m);
    components_DĪ“Ā²šā = zeros(eltype(U_jā),(deg+3)*m);
    DxĀ²Fā_eval = zeros(eltype(U_jā),deg+3); DxĪ“Fā_eval = zeros(eltype(U_jā),deg+3); 
    DĪ“Ā²Fā_eval = zeros(eltype(U_jā),deg+3);
    u = similar(component(U_jā,1));
    for j ā 1:m
        u[:] = coefficients(component(U_jā,j));
        for āā ā 1:deg+3
            DxĀ²Fā_eval[āā], DxĪ“Fā_eval[āā], DĪ“Ā²Fā_eval[āā] = šĀ²š( (j-1)/mi + (1/mi)*(nodes[āā]+1)/2, component(U_jā,j)(nodes[āā])[0] ,Ī“ ,Ī¦ ,para; Ī“ā)
        end
        components_DxĀ²šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DxĀ²Fā_eval,deg,degi,ppi);
        components_DxĪ“šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DxĪ“Fā_eval,deg,degi,ppi);
        components_DĪ“Ā²šā[1+(j-1)*(deg+3):j*(deg+3)] = cheb_interp(DĪ“Ā²Fā_eval,deg,degi,ppi);
    end
    šxĀ²šā = Sequence(Chebyshev(deg)^m,components_DxĀ²šā);
    šxĪ“šā = Sequence(Chebyshev(deg)^m,components_DxĪ“šā);
    šĪ“Ā²šā = Sequence(Chebyshev(deg)^m,components_DĪ“Ā²šā);
    return [šxĀ²šā šxĪ“šā šĪ“Ā²šā] 
end

# Check lag inclusion of numerical solution
function lag_range(u_tuple,Ī“,para,error;Ī“ā=1)
    lag = interpolate_lag(u_tuple,Ī“,para);
    u_jā,U_jā = u_tuple;    m = size(u_jā,2);
    ercheb = similar(lag(zero(eltype(U_jā))));  ercheb[:] .= para[4]*error;
    lag_eval = lag(Interval(-1,1)) - ercheb;
    range_lag = Interval(-1,1)*ones(Interval,m);
    for jā1:m
        maxv = sup(lag_eval[j]); minv = inf(lag_eval[j]);
        max_translate = (2/Ī“ā)*Interval(maxv)+1;
        min_translate = (2/Ī“ā)*Interval(minv)+1;
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
function bound_Y(A,G_vect,rā;do_round=0)
    if do_round==0
        Y = prenorm_vec(A*G_vect,rā);
    elseif do_round==1
        G! = interval.(Float64.(inf.(G_vect),RoundDown),Float64.(sup.(G_vect),RoundUp));
        Y = prenorm_vec(A*G!,rā);
    end
    return Y
end

# Yā bound
@views function bound_Yā(u_tuple,Ī“,Ī¦::Sequence{CartesianPower{Chebyshev}},para;Ī“ā=1,sdiv=10)   
    bstrap = 1;
    u_jā, U_jā = u_tuple;   k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = Interval(m);   Ī = 1/mi;
    bound_dįµF = zeros(eltype(U_jā),m);
    U = Array{Sequence}(undef,k+1);  U[1] = U_jā;
    U_eval = zeros(eltype(U_jā),k+1)
    scale_dĪ¦ = ((2/Ī“ā) .^ (0:k))*one(eltype(Ī¦));     #derivative scaling for Ī¦
    scale_dU = (2*m) .^ (0:k)*one(eltype(U_jā));     #derivative scaling for U
    for n=1:k+1-bstrap
        U[n+1] = Derivative(1)*U[n];
    end
    isd = interval(sdiv)
    for jā1:m
        for jj=1:sdiv
            tā±¼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tā±¼_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Ī¦_arg = Ī“*tā±¼ - para[3] - para[4]*component(U[1],j)(tā±¼_cheb)[0];
            Ī¦_arg = 1 .+ 2/Ī“ā*Ī¦_arg;    #rescale to correct domain [-1,1]
            for e=1:k+1
                U_eval[e] = component(U[e],j)(tā±¼_cheb)[0];
            end
            bound_dįµF[j] = interval( max( sup(bound_dįµF[j]), sup(abs(Fā(scale_dU.*U_eval[1:k+1],scale_dĪ¦.*coefficients(component(Ī¦,1:k+1)(Ī¦_arg)),para,Ī“;order=k)[k+1])) ) );
        end
    end
    Cā = 1/(Interval(factorial(k+1))*2^(2*(k+bstrap-1)));
    Yā = Ī^(k+bstrap)*Cā*norm(sup.(bound_dįµF),Inf);
    return Interval(sup(Yā))
end

@views function bound_Yā(u_tuple,Ī“,Ī¦,para;Ī“ā=1,sdiv=10)   
    bstrap = 1;
    u_jā, U_jā = u_tuple;   k = length(component(U_jā,1))-2;
    m = size(u_jā,2);   mi = Interval(m);   Ī = 1/mi;
    bound_dįµF = zeros(eltype(U_jā),m);
    U = Array{Sequence}(undef,k+1);  U[1] = U_jā;
    U_eval = zeros(eltype(U_jā),k+1)
    scale_dĪ¦ = ((2/Ī“ā) .^ (0:k))*one(eltype( Ī¦(interval(0))[1] ));     #derivative scaling for Ī¦
    scale_dU = (2*m) .^ (0:k)*one(eltype(U_jā));     #derivative scaling for U
    for n=1:k+1-bstrap
        U[n+1] = Derivative(1)*U[n];
    end
    isd = interval(sdiv)
    for jā1:m
        for jj=1:sdiv
            tā±¼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tā±¼_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Ī¦_arg = Ī“*tā±¼ - para[3] - para[4]*component(U[1],j)(tā±¼_cheb)[0];
            Ī¦_arg = 1 .+ 2/Ī“ā*Ī¦_arg;    #rescale to correct domain [-1,1]
            for e=1:k+1
                U_eval[e] = component(U[e],j)(tā±¼_cheb)[0];
            end
            bound_dįµF[j] = interval( max( sup(bound_dįµF[j]), sup(abs(Fā(scale_dU.*U_eval[1:k+1],scale_dĪ¦.*Ī¦(Ī¦_arg)[1:k+1],para,Ī“;order=k)[k+1])) ) );
        end
    end
    Cā = 1/(Interval(factorial(k+1))*2^(2*(k+bstrap-1)));
    Yā = Ī^(k+bstrap)*Cā*norm(sup.(bound_dįµF),Inf);
    return Interval(sup(Yā))
end

# Z0 bound
function bound_Z0(A,DG,rā;do_round=0)
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
    Z0 = prenorm_op(defect.coefficients,rā);
    return Z0
end

# Z1 bound
function bound_Z1(A,u_tuple,Ī“,s_jā,Ī¦,para,ppi,rā,r_ā; Ī“ā=1,sdiv=10)
    bstrap = 1;
    u_jā, U_jā = u_tuple;   k = size(u_jā,1)-2;    ki = Interval(k);
    m = size(u_jā,2);
    Ļ = zeros(eltype(u_jā),k+2,m);
    fac_p = one(eltype(u_jā))*factorial(bstrap);
    for āā1:k+2
        for jā1:m
            t = range(inf(s_jā[ā,j]) , sup(s_jā[1,j]) ;length=sdiv);
            DF_bound = zero(eltype(u_jā));
            for sd = 1:sdiv-1
                t_sd = @interval(t[sd])āŖ@interval(t[sd+1]);
                DF_bound = max( DF_bound, abs( šš(t_sd,component(U_jā,j)(-1+2*t_sd)[0], Ī“, Ī¦,para,1 ;Ī“ā )[1] ) )
            end
            Ļ[ā,j] = (abs(s_jā[ā,1] - 0)^bstrap)/fac_p * DF_bound;
        end
    end
    Ļ_vect = reshape(Ļ,(k+2)*(m),1);
    AĻ = [abs.(view(A,1:(k+2)*m,1:(k+2)*m))*Ļ_vect;0];
    return r_ā*prenorm_vec(AĻ,rā)
end

# Z2 bound
function bound_Z2(A,u_tuple,Ī“,s_jā,Ī¦::Sequence{CartesianPower{Chebyshev}},para,ppi,r,rā,r_ā;Ī“ā=1)
    bstrap = 1;
    u_jā, _ = u_tuple;   k = size(u_jā,1)-2;    m = size(u_jā,2); ki = Interval(k);
    Ī = 1 + (2\ppi)*log(ki+1);  mi = interval(m);
    Ļ±ā = Interval(0);   # Note: the hessian of the lag function vanishes FOR OUR EXAMPLE. This line of code is NOT GENERAL.
    Ļ±ā = zeros(eltype(u_jā),k+2,m); u_rad = r*(Ī+r_ā);
    Ī¦_bound = zeros(eltype(u_jā),3);    Ī“_rad = r*rā*Interval(-1,1);
    range_lag,_ = lag_range(u_tuple,Ī“+Ī“_rad,para,u_rad*Interval(-1,1); Ī“ā);   ## MOVE THIS COMPUTATION OUTSIDE OF THE Z2 BOUND.
    # Get inclusion for Ī¦
    for nā1:3
        for jā1:m
            Ī¦_bound[n] = Interval(max(sup(Ī¦_bound[n]),sup(component(Ī¦,n)(range_lag[j])[0])));
        end
    end
    bound_2derivative = zeros(eltype(u_jā),3,bstrap);
    for qā1:bstrap
        for j=1:m
            u_enc = Interval(-1,1)*(Ī*norm(u_jā[:,j],Inf) + u_rad);
            bound_2derivative[:,q] .= max.( bound_2derivative[:,q], DĀ²F(Interval(inf(s_jā[1,j]),sup(s_jā[end,j])),u_enc,Ī“+Ī“_rad,Interval(-1,1)*Ī¦_bound,para) );
        end
    end
    Hess_Ī = zeros(eltype(u_jā),2,2);
    for qā0:bstrap-1
        Hess_Ī[:] = abs.([bound_2derivative[1,q+1] bound_2derivative[2,q+1]; 
                    bound_2derivative[2,q+1] bound_2derivative[3,q+1] ]);
    end
    Īs_jā = similar(s_jā);
    for āā1:k+2
        Īs_jā[ā,:] = s_jā[ā,:] - s_jā[1,:]
    end
    int_multi = dot(Interval(-1,1)*[Ī+r_ā;rā],Hess_Ī*(Interval(-1,1)*[Ī+r_ā;rā]))/factorial(bstrap)
    int_bound = Īs_jā.^bstrap*int_multi;
    Ļ±ā = reshape(int_bound,m*(k+2));
    Ļ± = [Ļ±ā;Ļ±ā];
    Z2 = prenorm_vec(A*Ļ±*r,rā);
    return Z2
end

function bound_Z2(A,u_tuple,Ī“,s_jā,Ī¦,para,ppi,r,rā,r_ā;Ī“ā=1)
    bstrap = 1;
    u_jā, _ = u_tuple;   k = size(u_jā,1)-2;    m = size(u_jā,2); ki = Interval(k);
    Ī = 1 + (2\ppi)*log(ki+1);  mi = interval(m);
    Ļ±ā = Interval(0);   # Note: the hessian of the lag function vanishes FOR OUR EXAMPLE. This line of code is NOT GENERAL.
    Ļ±ā = zeros(eltype(u_jā),k+2,m); u_rad = r*(Ī+r_ā);
    Ī¦_bound = zeros(eltype(u_jā),3);    Ī“_rad = r*rā*Interval(-1,1);
    range_lag,_ = lag_range(u_tuple,Ī“+Ī“_rad,para,u_rad*Interval(-1,1); Ī“ā);   ## MOVE THIS COMPUTATION OUTSIDE OF THE Z2 BOUND.
    # Get inclusion for Ī¦
    for nā1:3
        for jā1:m
            Ī¦_bound[n] = Interval(max(sup(Ī¦_bound[n]),Ī¦(range_lag[j])[n]));
        end
    end
    bound_2derivative = zeros(eltype(u_jā),3,bstrap);
    for qā1:bstrap
        for j=1:m
            u_enc = Interval(-1,1)*(Ī*norm(u_jā[:,j],Inf) + u_rad);
            bound_2derivative[:,q] .= max.( bound_2derivative[:,q], DĀ²F(Interval(inf(s_jā[1,j]),sup(s_jā[end,j])),u_enc,Ī“+Ī“_rad,Interval(-1,1)*Ī¦_bound,para) );
        end
    end
    Hess_Ī = zeros(eltype(u_jā),2,2);
    for qā0:bstrap-1
        Hess_Ī[:] = abs.([bound_2derivative[1,q+1] bound_2derivative[2,q+1]; 
                    bound_2derivative[2,q+1] bound_2derivative[3,q+1] ]);
    end
    Īs_jā = similar(s_jā);
    for āā1:k+2
        Īs_jā[ā,:] = s_jā[ā,:] - s_jā[1,:]
    end
    int_multi = dot(Interval(-1,1)*[Ī+r_ā;rā],Hess_Ī*(Interval(-1,1)*[Ī+r_ā;rā]))/factorial(bstrap)
    int_bound = Īs_jā.^bstrap*int_multi;
    Ļ±ā = reshape(int_bound,m*(k+2));
    Ļ± = [Ļ±ā;Ļ±ā];
    Z2 = prenorm_vec(A*Ļ±*r,rā);
    return Z2
end

# Zā bound
function bound_Zā(u_tuple,Ī“,Ī¦::Sequence{CartesianPower{Chebyshev}},para,ppi,r,rā,r_ā; Ī“ā=1, sdiv=10)
    bstrap = 1;
    u_jā, U_jā = u_tuple;   k = size(u_jā,1)-2;    ki = Interval(k);
    Ī = 1 + (2\ppi)*log(ki+1);
    Cāā_star_1 = (1+Ī)*(ppi/4)^bstrap*factorial(k+1-bstrap)/Interval(factorial(k+1));
    Cāā_star_2 = 1/(Interval(factorial(bstrap))^2^bstrap)   # Note: always have kā„p and pā¤2.
    Cāā = Interval(min(sup(Cāā_star_1),sup(Cāā_star_2)));
    m = size(u_jā,2);    mi = Interval(m);
    Ī“_rad = r*rā*Interval(-1,1);    u_rad = r*(Ī+r_ā);
    u_perturb = u_rad*Interval(-1,1);
    bounds_DĪ = zeros(eltype(U_jā),m);
    scale_dĪ¦ = ((2/Ī“ā) .^ (0:1))*one(eltype(Ī¦));  
    isd = interval(sdiv);
    for jā1:m
        for jj=1:sdiv
            tā±¼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tā±¼_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Ī¦_arg = Ī“*tā±¼ - para[3] - para[4]*component(U_jā,j)(tā±¼_cheb)[0];
            Ī¦_arg = 1 .+ 2/Ī“ā*Ī¦_arg;    #rescale to correct domain [-1,1]
            U_eval = component(U_jā,j)(tā±¼_cheb)[0];
            šš1, šš2 = DF(tā±¼,U_eval + u_perturb,Ī“ + Ī“_rad,scale_dĪ¦.*coefficients(Ī¦(Ī¦_arg)),para);
            bounds_DĪ[j] = abs(šš1)*(Ī+r_ā) + abs(šš2)*(Ī+r_ā);
        end
    end
    Zā = Cāā*(1/mi)^bstrap*norm(bounds_DĪ,Inf);
    return Zā
end

function bound_Zā(u_tuple,Ī“,Ī¦,para,ppi,r,rā,r_ā; Ī“ā=1, sdiv=10)
    bstrap = 1;
    u_jā, U_jā = u_tuple;   k = size(u_jā,1)-2;    ki = Interval(k);
    Ī = 1 + (2\ppi)*log(ki+1);
    Cāā_star_1 = (1+Ī)*(ppi/4)^bstrap*factorial(k+1-bstrap)/Interval(factorial(k+1));
    Cāā_star_2 = 1/(Interval(factorial(bstrap))^2^bstrap) 
    Cāā = Interval(min(sup(Cāā_star_1),sup(Cāā_star_2)));
    m = size(u_jā,2);    mi = Interval(m);
    Ī“_rad = r*rā*Interval(-1,1);    u_rad = r*(Ī+r_ā);
    u_perturb = u_rad*Interval(-1,1);
    bounds_DĪ = zeros(eltype(U_jā),m);
    scale_dĪ¦ = ((2/Ī“ā) .^ (0:1))*one(eltype(Ī¦(interval(0))[1]));  
    isd = interval(sdiv);
    for jā1:m
        for jj=1:sdiv
            tā±¼ = interval(inf((j-1)/mi + (jj-1)/sdiv*(1/mi)),sup((j-1)/mi + (jj)/sdiv*(1/mi)))
            tā±¼_cheb = interval(inf(-1+2*(jj-1)/isd),sup(-1 + 2*jj/isd));
            Ī¦_arg = Ī“*tā±¼ - para[3] - para[4]*component(U_jā,j)(tā±¼_cheb)[0];
            Ī¦_arg = 1 .+ 2/Ī“ā*Ī¦_arg;    #rescale to correct domain [-1,1]
            U_eval = component(U_jā,j)(tā±¼_cheb)[0];
            šš1, šš2 = DF(tā±¼,U_eval + u_perturb,Ī“ + Ī“_rad,scale_dĪ¦.*Ī¦(Ī¦_arg)[1:2],para);
            bounds_DĪ[j] = abs(šš1)*(Ī+r_ā) + abs(šš2)*(Ī+r_ā);
        end
    end
    Zā = Cāā*(1/mi)^bstrap*norm(bounds_DĪ,Inf);
    return Zā
end

# radii polynomial
function radpol(G,A,DG,u,Ī“,Ī¦,para,s_jā,ppi,r_star,rā,r_ā;Ī“ā=1,rounding_YZ=0,debug=0,sdiv=10)
    if rounding_YZ==1
        u_jā! = interval.(Float64.(inf.(u[1]),RoundDown),Float64.(sup.(u[1]),RoundUp));
        U_jā! = Sequence(space(u[2]),interval.(Float64.(inf.(coefficients(u[2])),RoundDown),Float64.(sup.(coefficients(u[2])),RoundUp)));
        u! = (u_jā!,U_jā!);
        Ī“! = interval(Float64(inf(Ī“),RoundDown),Float64(sup(Ī“),RoundUp))
        Ī¦! = Sequence(space(Ī¦),interval.(Float64.(inf.(coefficients(Ī¦)),RoundDown),Float64.(sup.(coefficients(Ī¦)),RoundUp)));
        s_jā! = interval.(Float64.(inf.(s_jā),RoundDown),Float64.(sup.(s_jā),RoundUp));
        ppi! = interval(Float64(inf(ppi),RoundDown),Float64(sup(ppi),RoundUp))
        para! = interval.(Float64.(inf.(para),RoundDown),Float64.(sup.(para),RoundUp));
        Ī“ā! = interval(Float64(inf(Ī“ā),RoundDown),Float64(sup(Ī“ā),RoundUp))
        r_star! = interval(Float64(inf(r_star),RoundDown),Float64(sup(r_star),RoundUp))
        rā! = interval(Float64(inf(rā),RoundDown),Float64(sup(rā),RoundUp))
        r_ā! = interval(Float64(inf(r_ā),RoundDown),Float64(sup(r_ā),RoundUp))
    elseif rounding_YZ==0
        u! = u; Ī“! = Ī“; Ī¦! = Ī¦; s_jā! = s_jā;   ppi! = ppi; para! = para;   Ī“ā! = Ī“ā;   r_star! = r_star;   rā! = rā;   r_ā! = r_ā;
    end
    if debug==0
        Y0 = bound_Y(A,G,rā!;do_round=rounding_YZ);   sY = norm(sup.(Y0),Inf)
        println("   ā Yā bound computed: $sY.")
        Yā = bound_Yā(u!,Ī“!,Ī¦!,para!;Ī“ā,sdiv);   sYā = sup(Yā);
        println("   ā Yā bound computed: $sYā.")
        Z0 = bound_Z0(A,DG,rā!;do_round=rounding_YZ); sZ0 = norm(sup.(Z0),Inf);
        println("   ā Zā bound computed: $sZ0.")
        Z1 = bound_Z1(A,u!,Ī“!,s_jā!,component(Ī¦!,1:2),para!,ppi!,rā!,r_ā!;Ī“ā=Ī“ā!,sdiv);  sZ1 = norm(sup.(Z1),Inf);
        println("   ā Zā bound computed: $sZ1.")
        Z2 = bound_Z2(A,u!,Ī“!,s_jā!,component(Ī¦!,1:3),para!,ppi!,r_star!,rā!,r_ā!;Ī“ā=Ī“ā!);   sZ2 = norm(sup.(Z2),Inf);
        println("   ā Zā bound computed: $sZ2.")
        Zā = bound_Zā(u!,Ī“!,component(Ī¦!,1:2),para!,ppi!,r_star!,rā!,r_ā!;Ī“ā=Ī“ā!,sdiv);  sZā = sup(Zā);
        println("   ā Zā bound computed: $sZā.")
        poly_finite = [Y0 (Z2+Z1+Z0 .- 1)];
        poly_ā = [norm(Yā,Inf) Zā - r_ā];
        r = [- poly_finite[:,1]./poly_finite[:,2] ; -poly_ā[1]/poly_ā[2] ];
        if minimum(inf.(r))<0
            println("Existence interval is empty; a radius is negative.")
            ie = ā
            C0_err = Inf;   Ī“_err = Inf;
            # @infiltrate
        elseif maximum(nextfloat.(sup.(r)))>r_star
            println("Existence interval is empty; a radius is larger than r_star.")
            ie = ā
            C0_err = Inf;   Ī“_err = Inf;
            # @infiltrate
        else
            ie = Interval(maximum(nextfloat.(sup.(r))),inf(r_star))
            k = size(u[1],1)-2;    ki = Interval(k);
            Ī = 1 + (2\ppi)*log(ki+1);
            C0_err = inf(ie)*(Ī+r_ā); Ī“_err = inf(ie)*rā;   sC0 = sup(C0_err);  sĪ“ = sup(Ī“_err);
            println("   ā Cā° enclosure: $sC0.")
            println("   ā Ī“ enlosure: $sĪ“.")
        end
        return Y0, Yā, Z0, Z1, Z2, Zā, poly_finite, poly_ā, r, ie, C0_err, Ī“_err
    else
        # @infiltrate
        return []
    end
end

function radpol(G,A,DG,u,Ī“,Ī¦_Cheb,Ī¦_function,para,s_jā,ppi,r_star,rā,r_ā;Ī“ā=1,rounding_YZ=0,debug=0,sdiv=10)
    if rounding_YZ==1
        u_jā! = interval.(Float64.(inf.(u[1]),RoundDown),Float64.(sup.(u[1]),RoundUp));
        U_jā! = Sequence(space(u[2]),interval.(Float64.(inf.(coefficients(u[2])),RoundDown),Float64.(sup.(coefficients(u[2])),RoundUp)));
        u! = (u_jā!,U_jā!);
        Ī“! = interval(Float64(inf(Ī“),RoundDown),Float64(sup(Ī“),RoundUp))
        Ī¦_Cheb! = Sequence(space(Ī¦_Cheb),interval.(Float64.(inf.(coefficients(Ī¦_Cheb)),RoundDown),Float64.(sup.(coefficients(Ī¦_Cheb)),RoundUp)));
        s_jā! = interval.(Float64.(inf.(s_jā),RoundDown),Float64.(sup.(s_jā),RoundUp));
        ppi! = interval(Float64(inf(ppi),RoundDown),Float64(sup(ppi),RoundUp))
        para! = interval.(Float64.(inf.(para),RoundDown),Float64.(sup.(para),RoundUp));
        Ī“ā! = interval(Float64(inf(Ī“ā),RoundDown),Float64(sup(Ī“ā),RoundUp))
        r_star! = interval(Float64(inf(r_star),RoundDown),Float64(sup(r_star),RoundUp))
        rā! = interval(Float64(inf(rā),RoundDown),Float64(sup(rā),RoundUp))
        r_ā! = interval(Float64(inf(r_ā),RoundDown),Float64(sup(r_ā),RoundUp))
        Ī¦_function! = t-> interval.(Float64.(inf.(Ī¦_function(t)),RoundDown),Float64.(sup.(Ī¦_function(t)),RoundUp));
    elseif rounding_YZ==0
        u! = u; Ī“! = Ī“; Ī¦_Cheb! = Ī¦_Cheb; s_jā! = s_jā;   ppi! = ppi; para! = para;   Ī“ā! = Ī“ā; r_star! = r_star;   rā! = rā;   r_ā! = r_ā; Ī¦_function! = Ī¦_function;
    end
    if debug==0
        Y0 = bound_Y(A,G,rā!;do_round=rounding_YZ);   sY = norm(sup.(Y0),Inf)
        println("   ā Yā bound computed: $sY.")
        Yā = bound_Yā(u!,Ī“!,Ī¦_function!,para!;Ī“ā,sdiv);   sYā = sup(Yā);
        println("   ā Yā bound computed: $sYā.")
        Z0 = bound_Z0(A,DG,rā;do_round=rounding_YZ); sZ0 = norm(sup.(Z0),Inf);
        println("   ā Zā bound computed: $sZ0.")
        Z1 = bound_Z1(A,u!,Ī“!,s_jā!,component(Ī¦_Cheb!,1:2),para!,ppi!,rā!,r_ā!;Ī“ā=Ī“ā!,sdiv);  sZ1 = norm(sup.(Z1),Inf);
        println("   ā Zā bound computed: $sZ1.")
        Z2 = bound_Z2(A,u!,Ī“!,s_jā!,Ī¦_function!,para!,ppi!,r_star!,rā!,r_ā!;Ī“ā=Ī“ā!);   sZ2 = norm(sup.(Z2),Inf);
        println("   ā Zā bound computed: $sZ2.")
        Zā = bound_Zā(u!,Ī“!,Ī¦_function!,para!,ppi!,r_star!,rā!,r_ā!;Ī“ā=Ī“ā!,sdiv);  sZā = sup(Zā);
        println("   ā Zā bound computed: $sZā.")
        poly_finite = [Y0 (Z2+Z1+Z0 .- 1)];
        poly_ā = [norm(Yā,Inf) Zā - r_ā];
        r = [- poly_finite[:,1]./poly_finite[:,2] ; -poly_ā[1]/poly_ā[2] ];
        if minimum(inf.(r))<0
            println("Existence interval is empty; a radius is negative.")
            ie = ā
            C0_err = Inf;   Ī“_err = Inf;
            # @infiltrate
        elseif maximum(nextfloat.(sup.(r)))>r_star
            println("Existence interval is empty; a radius is larger than r_star.")
            ie = ā
            C0_err = Inf;   Ī“_err = Inf;
            # infiltrate
        else
            ie = Interval(maximum(nextfloat.(sup.(r))),inf(r_star))
            k = size(u[1],1)-2;    ki = Interval(k);
            Ī = 1 + (2\ppi)*log(ki+1);
            C0_err = inf(ie)*(Ī+r_ā); Ī“_err = inf(ie)*rā;   sC0 = sup(C0_err);  sĪ“ = sup(Ī“_err);
            println("   ā Cā° enclosure: $sC0.")
            println("   ā Ī“ enlosure: $sĪ“.")
        end
        return Y0, Yā, Z0, Z1, Z2, Zā, poly_finite, poly_ā, r, ie, C0_err, Ī“_err
    else
        # @infiltrate
        return []
    end
end

function modify_candidate_zero(iu,Ī“,para,Ī“ā;desired_tolerance=1e-20)
    Ī± = para[3];    c = para[4]
    u_0 = iu[1][1,1];    u_1 = iu[1][end,end]
    scale = sup.( [Ī“ -c*u_1 ; 0 -c*u_0]\[-desired_tolerance - (Ī“ - Ī± - c*u_1) ; -Ī“ā + desired_tolerance - (-Ī± - c*u_0)] )
    while sup((1+scale[1])*Ī“ - Ī± - (1+scale[2])*c*u_1)>0 && inf((1/Ī“ā)*(-Ī± - (1+scale[2])*c*u_0))<-1
        scale = scale*1.01;
    end
    iu_mod = ( interval(1+scale[2])*iu[1], Sequence(space(iu[2]),interval(1+scale[2])*coefficients(iu[2])) );
    Ī“_mod = (1+scale[2])*Ī“;
    return iu_mod,Ī“_mod
end

function modify_candidate_zero_rightside_only(iu,Ī“,para;desired_tolerance=1e-20)
    Ī± = para[3];    c = para[4]
    u_1 = iu[1][end,end]
    if sup(Ī“ - Ī± - c*u_1)>0
        Ļµā = (-interval(desired_tolerance) - (Ī“ - Ī± - c*u_1))/Ī“;
        while sup((1+Ļµā)*Ī“ - Ī± -c*u_1)>0
            Ļµā = Ļµā*1.01;
        end
        Ī“_mod = (1+Ļµā)*Ī“
    else
        Ī“_mod = Ī“
    end
    return iu,Ī“_mod
end

function check_lag_monotonicity(Ī¦::Sequence{CartesianPower{Chebyshev}},U,s_jā,Ī“,r_C0)
    # Note 1, this code is example-specific.
    Ī³,Īŗ,Ī±,c = para
    t = LinRange(-1,1,sdiv)
    d_lag = Inf
    du = interval.(zero(eltype(U)))
    u = interval.(zero(eltype(U)))
    Ī¦_eval = interval.(zero(eltype(U)))
    dĪ¦_eval = interval.(zero(eltype(U)))
    for n=1:num_components(U)
        u = component(U,n)(interval(-1,1))[0] + r_C0*interval(-1,1)
        t = interval(inf(s_jā[1,n]),sup(s_jā[end,n]))
        Ī¦_eval = component(Ī¦,1)(t*Ī“ - Ī± - c*u)[0]
        dĪ¦_eval = component(Ī¦,2)(t*Ī“ - Ī± - c*u)[0]
        du = Ī“*(-Ī³*u - Īŗ*Ī¦_eval)
        d_lag = min(d_lag, inf(Ī“ - c*du) )
    end
    check = d_lag>0
    println("   ā The lag is monotone: $check.")
    return d_lag>0
end

function check_lag_monotonicity(para,Ī¦,U,s_jā,Ī“,r_C0)
    # Note 1, this code is example-specific.
    Ī³,Īŗ,Ī±,c = para
    d_lag = Inf
    du = interval.(zero(eltype(U)))
    u = interval.(zero(eltype(U)))
    Ī¦_eval = interval.(zero(eltype(U)))
    dĪ¦_eval = interval.(zero(eltype(U)))
    for n=1:num_components(U)
        u = component(U,n)(interval(-1,1))[0] + r_C0*interval(-1,1)
        t = interval(inf(s_jā[1,n]),sup(s_jā[end,n]))
        Ī¦_eval = Ī¦(t*Ī“ - Ī± - c*u)[1]
        dĪ¦_eval = Ī¦(t*Ī“ - Ī± - c*u)[2]
        du = Ī“*(-Ī³*u - Īŗ*Ī¦_eval)
        d_lag = min(d_lag, inf(Ī“ - c*du) )
    end
    check = d_lag>0
    println("   ā The lag is monotone: $check.")
    return d_lag>0
end

function proof_section_radpol_interp_monotonicity(G_n,A_n,DG_n,iu_n,iĪ“_n,Ī¦,Ī¦_function,ipara,s_jā_n,ppi,r_star,r_0,r_ā,Ī“_val_previous,kĪ¦,kĪ¦_low,š¤,rounding_YZ)
    println(": Starting evaluation of the bounds, radii polynomials. :")
    _, _, _, _, _, _, _, _, _, _, C0_err, Ī“_err = radpol(G_n,A_n,DG_n,iu_n,iĪ“_n,Ī¦,Ī¦_function,ipara,s_jā_n,ppi,r_star,r_0,r_ā;Ī“ā=Ī“_val_previous,rounding_YZ);
    println(": Starting high-accuracy interpolation. :")
    _,_,Ī¦_interp_high_accuracy = interpolate_solution_tight_D012(iu_n[2],iĪ“_n,s_jā_n,Ī¦_function,C0_err,Ī“_err,kĪ¦,interval(kĪ¦),š¤,ipara,ppi;Ī“ā=Ī“_val_previous,sdiv=10,check_large_coeffs=1)
    println(": Starting low-accuracy interpolation. :")
    _,_,Ī¦_interp_low_accuracy = interpolate_solution_tight_D012(iu_n[2],iĪ“_n,s_jā_n,Ī¦_function,C0_err,Ī“_err,kĪ¦_low,interval(kĪ¦_low),š¤,ipara,ppi;Ī“ā=Ī“_val_previous,sdiv=10,max_N=kĪ¦_low,check_large_coeffs=0)
    println(": Building hybrid enclosure. :")
    Ī¦_function_next = t -> evaluation_hybrid_enclosure(t,convert_Chebyshev_to_callable_function(component(Ī¦_interp_high_accuracy,1)),
                            iĪ“_n + interval(-1,1)*Ī“_err,Ī¦_function,ipara,all_derivatives(iu_n[1],ipara,Ī¦_function,iĪ“_n + interval(-1,1)*Ī“_err;
                            Ī“ā=Ī“_val_previous);n_derivatives=š¤,Ī“ā=Ī“_val_previous)
    println(": Verifying monotonicity of time lag. :")
    check_lag_monotonicity(ipara,Ī¦_function,iu_n[2],s_jā_n,iĪ“_n + interval(-1,1)*Ī“_err,C0_err)
    return C0_err, Ī“_err, Ī¦_interp_high_accuracy, Ī¦_interp_low_accuracy, Ī¦_function_next
end

function proof_section_radpol_interp_monotonicity_nohybrid(G_n,A_n,DG_n,iu_n,iĪ“_n,Ī¦,Ī¦_function,ipara,s_jā_n,ppi,r_star,r_0,r_ā,Ī“_val_previous,kĪ¦,kĪ¦_low,š¤,rounding_YZ)
    println(": Starting evaluation of the bounds, radii polynomials. :")
    _, _, _, _, _, _, _, _, _, _, C0_err, Ī“_err = radpol(G_n,A_n,DG_n,iu_n,iĪ“_n,Ī¦,Ī¦_function,ipara,s_jā_n,ppi,r_star,r_0,r_ā;Ī“ā=Ī“_val_previous,rounding_YZ);
    println(": Starting high-accuracy interpolation. :")
    _,_,Ī¦_interp_high_accuracy = interpolate_solution_tight_D012(iu_n[2],iĪ“_n,s_jā_n,Ī¦_function,C0_err,Ī“_err,kĪ¦,interval(kĪ¦),š¤,ipara,ppi;Ī“ā=Ī“_val_previous,sdiv=10,check_large_coeffs=1)
    println(": Starting low-accuracy interpolation. :")
    _,_,Ī¦_interp_low_accuracy = interpolate_solution_tight_D012(iu_n[2],iĪ“_n,s_jā_n,Ī¦_function,C0_err,Ī“_err,kĪ¦_low,interval(kĪ¦_low),š¤,ipara,ppi;Ī“ā=Ī“_val_previous,sdiv=10,max_N=kĪ¦_low,check_large_coeffs=0)
    println(": Verifying monotonicity of time lag. :")
    check_lag_monotonicity(ipara,Ī¦_function,iu_n[2],s_jā_n,iĪ“_n + interval(-1,1)*Ī“_err,C0_err)
    return C0_err, Ī“_err, Ī¦_interp_high_accuracy, Ī¦_interp_low_accuracy
end
