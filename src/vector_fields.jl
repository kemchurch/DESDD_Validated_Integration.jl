function τ(t,u,p)
    return p[3] + p[4]*u
end

function τ(u,p)
    return p[3] + p[4]*u
end

function D₂τ(t,u,p)
    return p[4]
end

function Dx(x,ϕ,p,δ)
    γ,κ,α,c = p;
    F = δ*(-γ*x - κ*ϕ);
    return F
end

function D²x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx = X;
    ϕ,dϕ = κ*Φ;
    h1 = δ - c*dx;
    F = δ*(-γ*dx - dϕ*h1) ;
    return F
end

function D³x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x = X;
    ϕ,dϕ,d²ϕ = κ*Φ; 
    h1 = δ - c*dx;  h2 = -c*d²x;
    F = δ*(-γ*d²x - d²ϕ*h1^2 - dϕ*h2);
    return F
end

function D⁴x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x,d³x = X;
    ϕ,dϕ,d²ϕ,d³ϕ = κ*Φ;
    h1 = δ-c*dx;    h2 = -c*d²x;    h3 = -c*d³x;
    F = δ*(-γ*d³x - d³ϕ*h1^3 - dϕ*h3 - 3*d²ϕ*h1*h2);
    return F
end

function D⁵x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x,d³x,d⁴x = X;
    ϕ,dϕ,d²ϕ,d³ϕ,d⁴ϕ = κ*Φ;
    h1 = δ-c*dx;    h2 = -c*d²x;    h3 = -c*d³x;    h4 = -c*d⁴x;
    F = δ*(-γ*d⁴x - 3*d²ϕ*h2^2 - d⁴ϕ*h1^4 - dϕ*h4 - 4*d²ϕ*h1*h3 - 6*d³ϕ*h1^2*h2);
    return F
end

function D⁶x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x,d³x,d⁴x,d⁵x = X;
    ϕ,dϕ,d²ϕ,d³ϕ,d⁴ϕ,d⁵ϕ = κ*Φ;
    h1 = δ-c*dx;    h2 = -c*d²x;    h3 = -c*d³x;    h4 = -c*d⁴x;    h5 = -c*d⁵x;
    F = δ*( -γ*d⁵x - dϕ*h5 - d⁵ϕ*h1^5 - 10*d³ϕ*h1^2*h3 - 15*d³ϕ*h1*h2 - 10*d⁴ϕ*h1*h2 - 5*d²ϕ*h1*h4 - 10*d²ϕ*h2*h3 );
    return F
end

function D⁷x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x,d³x,d⁴x,d⁵x,d⁶x = X;
    ϕ,dϕ,d²ϕ,d³ϕ,d⁴ϕ,d⁵ϕ,d⁶ϕ = κ*Φ;
    h1 = δ-c*dx;    h2 = -c*d²x;    h3 = -c*d³x;    h4 = -c*d⁴x;    h5 = -c*d⁵x;    h6 = -c*d⁶x
    F = δ*( -γ*d⁶x - 10*d²ϕ*h3^2 - d⁶ϕ*h1^6 - dϕ*h6 - 15*d³ϕ*h2^3 - 15*d⁵ϕ*h1^4*h2 - 45*d⁴ϕ*h1^2*h2^2 - 20*d⁴ϕ*h1^3*h3 - 6*d²ϕ*h1*h5 - 15*d²ϕ*h2*h4 
            -15*d³ϕ*h1^2*h4 - 60*d³ϕ*h1*h2*h3 );
    return F
end

function D⁸x(X,Φ,p,δ)
    γ,κ,α,c = p;
    x,dx,d²x,d³x,d⁴x,d⁵x,d⁶x,d⁷x = X;
    ϕ,dϕ,d²ϕ,d³ϕ,d⁴ϕ,d⁵ϕ,d⁶ϕ,d⁷ϕ = κ*Φ;
    h1 = δ-c*dx;    h2 = -c*d²x;    h3 = -c*d³x;    h4 = -c*d⁴x;    h5 = -c*d⁵x;    h6 = -c*d⁶x;    h7 = -c*d⁷x;
    F = δ*( -γ*d⁶x - dϕ*h7 - d⁷ϕ*h1^7 - 21*d³ϕ*h1^2*h5 - 70*d³ϕ*h1*h3^2 - 21*d⁶ϕ*h1^5*h2 - 21*d⁶ϕ*h1^5*h2 - 105*d⁵ϕ*h1^3*h2^2 - 105*d³ϕ*h2^2*h3
            -35*d⁵ϕ*h1^4*h3 - 35*d⁴ϕ*h1^3*h4 - 7*d²ϕ*h1*h6 - 21*d²ϕ*h2*h5 - 35*d²ϕ*h3*h4 - 105*d⁴ϕ*h1*h2^3 - 105*d²ϕ*h1*h2*h4 - 210*d⁴ϕ*h1^2*h2*h3 );
    return F
end

function compute_max_degree(deg_X,deg_Φ)
    Deg_Φ = zeros(Float64,6);   Deg_Φ[1] = deg_Φ;
    Deg_Φ⁰ = max(deg_X,1);
    for n ∈ 2:6
        if Deg_Φ[n-1]==0 || Deg_Φ[n-1]==-Inf
            Deg_Φ[n] = -Inf;
            else
            Deg_Φ[n] = max(Deg_Φ[n-1]-1,0);
        end
    end
    deg_xdot1 = max(deg_X[1], Deg_Φ[1]*Deg_Φ⁰);
    deg_xdot2 = max(deg_xdot1, Deg_Φ[2]*Deg_Φ⁰+deg_xdot1);
    deg_xdot3 = max(deg_xdot2, Deg_Φ[3]*Deg_Φ⁰+deg_xdot1*2, Deg_Φ[2]*Deg_Φ⁰ + deg_xdot2);
    deg_xdot4 = max(deg_xdot3, Deg_Φ[4]*Deg_Φ⁰+deg_xdot1*3, Deg_Φ[3]*Deg_Φ⁰ + deg_xdot2 + deg_xdot1, Deg_Φ[2]*Deg_Φ⁰ + deg_xdot3);
    deg_xdot5 = max(deg_xdot4, Deg_Φ[5]*Deg_Φ⁰+deg_xdot1*4, Deg_Φ[4]*Deg_Φ⁰ + deg_xdot2 + deg_xdot1*2, Deg_Φ[3]*Deg_Φ⁰+ deg_xdot3 + deg_xdot1, Deg_Φ[3]*Deg_Φ⁰+deg_xdot2*2, Deg_Φ[3]*Deg_Φ⁰ + deg_xdot1 + deg_xdot3, Deg_Φ[2]*Deg_Φ⁰ + deg_xdot4);
    deg_xdot6 = max(deg_xdot5, Deg_Φ[5]*Deg_Φ⁰+deg_xdot1*5, deg_xdot4*Deg_Φ⁰ + deg_xdot2 + 4*deg_xdot1, deg_xdot4*Deg_Φ⁰ + deg_xdot2 + 3*deg_xdot1, deg_xdot3*Deg_Φ⁰ + deg_xdot3 + 2*deg_xdot1, deg_xdot2*Deg_Φ⁰ + deg_xdot1 + deg_xdot4, deg_xdot3*Deg_Φ⁰ + deg_xdot1 + 2*deg_xdot2, deg_xdot2*Deg_Φ⁰ + deg_xdot2 + deg_xdot2, deg_xdot1*Deg_Φ⁰ + deg_xdot5);
    return Int64.([deg_xdot1,deg_xdot2,deg_xdot3,deg_xdot4,deg_xdot5,deg_xdot6])
end

function compute_max_degree_DF(deg_Fₚ,deg_X,deg_Φ)
    Deg_Φ = zeros(Float64,6);   Deg_Φ[1] = deg_Φ;
    Deg_Φ⁰ = max(deg_X,1);
    for n ∈ 2:6
        if Deg_Φ[n-1]==0 || Deg_Φ[n-1]==-Inf
            Deg_Φ[n] = -Inf;
            else
            Deg_Φ[n] = max(Deg_Φ[n-1]-1,0);
        end
    end
    degr = deg_Fₚ + Deg_Φ[2]*Deg_Φ⁰;
    return Int64.(degr)
end

function compute_max_degree_D²F(deg_Fₚ,deg_X,deg_Φ)
    Deg_Φ = zeros(Float64,5);   Deg_Φ[1] = deg_Φ;
    Deg_Φ⁰ = max(deg_X,1);
    for n ∈ 2:6
        if Deg_Φ[n-1]==0 || Deg_Φ[n-1]==-Inf
            Deg_Φ[n] = -Inf;
            else
            Deg_Φ[n] = max(Deg_Φ[n-1]-1,0);
        end
    end
    degr = deg_Fₚ + 2*Deg_Φ[2]*Deg_Φ⁰;
    return Int64.(degr)
end

function F(x,Φ,p,δ)   # vector field evaluation
    dx  = Dx(x,Φ[1],p,δ);
    return dx
end

function F₈(x,Φ,p,δ;order=8)   # Time derivatives 0…7 of vector field.
    dx  = Dx(x[1],Φ[1],p,δ);
    if order==0
        return [dx]
    end
    d²x = D²x((x[1],x[2]),Φ[1:2],p,δ);
    if order==1
        return [dx;d²x]
    end
    d³x = D³x((x[1],x[2],x[3]),Φ[1:3],p,δ);
    if order==2
        return [dx;d²x;d³x]
    end
    d⁴x = D⁴x((x[1],x[1],x[3],x[4]),Φ[1:4],p,δ);
    if order==3
        return [dx;d²x;d³x;d⁴x]
    end
    d⁵x = D⁵x((x[1],x[1],x[3],x[4],x[5]),Φ[1:5],p,δ);
    if order==4
        return [dx;d²x;d³x;d⁴x;d⁵x]
    end
    d⁶x = D⁶x((x[1],x[1],x[3],x[4],x[5],x[6]),Φ[1:6],p,δ);
    if order==5
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x]
    end
    d⁷x = D⁷x((x[1],x[1],x[3],x[4],x[5],x[6],x[7]),Φ[1:7],p,δ);
    if order==6
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x;d⁷x]
    end
    d⁸x = D⁸x((x[1],x[1],x[3],x[4],x[5],x[6],x[7],x[8]),Φ,p,δ);
    if order==7
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x;d⁷x;d⁸x]
    end
    error("Error! only order ≤7 is supported!")
end

function F₈_bootstrap(x,Φ,p,δ;order=8)
    dx  = Dx(x,Φ[1],p,δ);
    if order==1
        return [dx]
    end
    d²x = D²x((x,dx),Φ[1:2],p,δ);
    if order==2
        return [dx;d²x]
    end
    d³x = D³x((x,dx,d²x),Φ[1:3],p,δ);
    if order==3
        return [dx;d²x;d³x]
    end
    d⁴x = D⁴x((x,dx,d²x,d³x),Φ[1:4],p,δ);
    if order==4
        return [dx;d²x;d³x;d⁴x]
    end
    d⁵x = D⁵x((x,dx,d²x,d³x,d⁴x),Φ[1:5],p,δ);
    if order==5
        return [dx;d²x;d³x;d⁴x;d⁵x]
    end
    d⁶x = D⁶x((x,dx,d²x,d³x,d⁴x,d⁵x),Φ[1:6],p,δ);
    if order==6
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x]
    end
    d⁷x = D⁷x((x,dx,d²x,d³x,d⁴x,d⁵x,d⁶x),Φ[1:7],p,δ);
    if order==7
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x;d⁷x]
    end
    d⁸x = D⁸x((x,dx,d²x,d³x,d⁴x,d⁵x,d⁶x,d⁷x),Φ,p,δ);
    if order==8
        return [dx;d²x;d³x;d⁴x;d⁵x;d⁶x;d⁷x;d⁸x]
    end
    error("Error! only order ≤8 is supported!")
end

function DF(t,x,δ,Φ,p)  # Derivatives of bootstrapped vector fields w.r.t. x, δ
    γ,κ,α,c = p;
    ϕ,dϕ = Φ[1:2];
        D₁F = (-γ - κ*dϕ*(-c))*δ;
        D₂F = -γ*x - κ*ϕ + t*(-κ*dϕ*δ);   
    return D₁F, D₂F
end

function D²F(t,x,δ,Φ,p);
    γ,κ,α,c = p;
    ϕ,dϕ,d²ϕ = Φ[1:3];
    Dₓ²F = (κ*c^2*d²ϕ)*δ
    D_δ²F = t^2*κ*d²ϕ;
    DₓD_δF = -γ + κ*c*dϕ - t*κ*c*δ*d²ϕ;
    return Dₓ²F, DₓD_δF, D_δ²F
end

function 𝐅(t,u,δ,Φ,p; δ₋= 1, check_edges=0)
    # Evaluation of bootstrapped vector fields, delta-scale, for specified argument t∈[0,1], u(t), given initial data Φ and parameter data p.
    Φ_arg = t*δ .- (p[3] + p[4]*u);     # argument of Φ in [-δ₋,0], interval'd.
    Φ_arg = 1 .+ 2/δ₋*Φ_arg;                      # translate to [-1,1]
    if check_edges == 0
        eval = Φ(Φ_arg)[0]
    elseif check_edges == 1                 
        if sup(interval.(Φ_arg))<-1
            DΦ = Derivative(1)*Φ;
            eval = Φ(-1)[0] + (DΦ(-1)[0])*(Φ_arg + 1);
        elseif inf(interval.(Φ_arg))<-1 && sup(interval.(Φ_arg))>=-1
            DΦ = Derivative(1)*Φ;
            split_left = interval(inf(interval.(Φ_arg)), -1);   split_right = interval(-1,sup(interval.(Φ_arg)))
            eval = union( Φ(-1)[0] + (DΦ(-1)[0])*(split_left + 1) ,  Φ(split_right)[0]  )
        elseif inf(interval(Φ_arg))>1
            DΦ = Derivative(1)*Φ;
            eval = Φ(1)[0] + (DΦ(1)[0])*(Φ_arg - 1);
        elseif inf(interval(Φ_arg))<=1 && sup(interval(Φ_arg))>1
            DΦ = Derivative(1)*Φ;
            split_left = interval(inf(interval.(Φ_arg)), 1);   split_right = interval(1,sup(interval.(Φ_arg)))
            eval = union( Φ(split_left)[0], Φ(1)[0] + (DΦ(-1)[0])*(split_right - 1) )
        else
            eval = Φ(Φ_arg)[0]
        end
    end
    return 𝐅 = (F(u,eval,p,δ));  # Compute and δ-scale
end

function 𝐅₈(t,u,δ,Φ::Sequence{CartesianPower{Chebyshev}},para,ind;  δ₋= 1, U_error=0, order=8, components_Φ=9)
    # Evaluation of bootstrapped vector fields, delta-scale, for specified argument t∈[0,1], u(t), given initial data Φ and parameter data p.
    Φ_arg = t*δ .- (para[3] .+ para[4]*u);         # argument of Φ in [-δ₋,0]
    Φ_arg = 1 .+ 2/δ₋*Φ_arg;                   # scale argument to [-1,1]
    scale = ((2/δ₋) .^ (0:components_Φ-1))*one(eltype(u)); 
    # compute derivative scale factor due to domain shift [-δ₋,0] to [-1,1]
    return F₈_bootstrap(u .+ U_error*interval(-1,1),scale.*Φ(Φ_arg),para,δ; order)[ind]; # Compute and δ-scale
end

function 𝐅₈(t,u,δ,Φ_function,para;  δ₋= 1, U_error=0, components_Φ=9)
    # Evaluation of bootstrapped vector fields, delta-scale, for specified argument t∈[0,1], u(t), given initial data Φ and parameter data p.
    # Note, in this scope, Φ_function is a function t↦[Φ₀(t);…;Φ₇(t)].
    Φ_arg = t*δ .- (para[3] .+ para[4]*u);         # argument of Φ in [-δ₋,0]
    Φ_arg = 1 .+ 2/δ₋*Φ_arg;                   # scale argument to [-1,1]
    scale = ((2/δ₋) .^ (0:components_Φ-1))*one(eltype(u));    # compute derivative scale factor due to domain shift [-δ₋,0] to [-1,1]
    return F₈_bootstrap(u .+ U_error*interval(-1,1),scale.*Φ_function(Φ_arg),para,δ); # Compute and δ-scale
end