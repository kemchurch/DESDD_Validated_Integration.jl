# Functions for interpolation and conversion

function cheb_nodes!(N::Integer)
    # Compute the Chebyshev nodes, rigorous enclosure.
    ppi = interval(pi); Ni = interval(N);
    req_nodes = cos.(ppi*(BigFloat.(0:1:N))/((Ni)));
    return req_nodes
end

function cheb_nodes(N::Integer,Ni::T where T<:Real,ppi::T where T<:Real)
    # Compute the Chebyshev nodes, rigorous enclosure, type stable provided N is an integer, Ni and ppi have the same type.
    req_nodes = cos.(ppi*((0:1:N))/((Ni)));
    return req_nodes
end

function cheb!(f::Any,N::Integer)
    # Convert generic function f to its associated Chebyshev coefficients. Slow due to type Any. Not used by the proofs. Function name is not correct Julian style.
    ppi = @interval pi; Ni = @interval N;
    req_nodes = cheb_nodes(N,Ni,ppi);
    f₀ = f.(req_nodes);
    a = interval.(zeros(BigFloat,N+1));
    
    e₁ = exp.(im*ppi*(0:1:N)/Ni);
    e₂ = exp.(-im*ppi*(0:1:N)/Ni);
    for n in 0:N
        T = real.(e₁.^n + e₂.^n)/2;
        T[1] = T[1]/2;  T[end] = T[end]/2;
        aval = 2*(T'*f₀)/Ni;
        a[n+1] = aval[1];
    end
    a[2:end] = a[2:end]/2;
    return a
end

function cheb_interp!(f,N::Integer)
    # Convert a generic function f to a Chebyshev interpolant. Outputs as a RadiiPolynomial object. Slow. No longer used by the proofs.
    a = cheb!(f,N);
    cf = Sequence(Chebyshev(N),a);
    return cf
end

function cheb_interp(f₀,N::Integer,Ni::T where T<:Real,ppi::T where T<:Real)
    # Type stable chebyshev interpolation, provided inputs are formatted properly. See cheb_nodes.
    a = zeros(eltype(f₀),N+1);
    e₁ = exp.(im*ppi*(0:1:N)/Ni);
    e₂ = exp.(-im*ppi*(0:1:N)/Ni);
    T = similar(a);
    for n in 0:N
        T[:] = real.(e₁.^n + e₂.^n)/2;
        T[1] = T[1]/2;  T[end] = T[end]/2;
        if n==N
            T[:] = T[:]/2;
        end
        a[n+1] = (T'*f₀)/Ni;
    end
    return a
end

@views function convert_matrix_to_interp(u,ppi::T where T<:Real)
    m = size(u,2);
    k = size(u,1)-1;    ki = one(eltype(u))*k;
    nodes = reverse(cheb_nodes(k,ki,ppi));
    ChebIdentity = Sequence(Chebyshev(k)^(k+1),one(eltype(u))*reshape(I[1:k+1,1:k+1],(k+1)^2));
    ChebMatrix_VectorOfVectors = coefficients.(ChebIdentity.(nodes));
    ChebMatrix = zeros(eltype(u),k+1,k+1);
    for j ∈ 1:k+1
        for i ∈ 1:k+1
            ChebMatrix[i,j] = ChebMatrix_VectorOfVectors[i][j];
        end
    end
    ChebMatrix_inv = inv(ChebMatrix);
    U_jℓ_entries = reshape(ChebMatrix_inv*u,(k+1)*m)
    U_jℓ = Sequence(Chebyshev(k)^m,U_jℓ_entries);
    return U_jℓ
end

function Cheb_node_to_spectral_matrix(u,nodes_input)
    one! = one(eltype(u))*one(eltype(nodes_input));
    k = size(u,1)-1;
    ChebIdentity = Sequence(Chebyshev(k)^(k+1),one!*reshape(I[1:k+1,1:k+1],(k+1)^2));
    ChebMatrix_VectorOfVectors = coefficients.(ChebIdentity.(nodes_input));
    ChebMatrix = zeros(eltype(u),k+1,k+1)*one!;
    for j ∈ 1:k+1
        for i ∈ 1:k+1
            ChebMatrix[i,j] = ChebMatrix_VectorOfVectors[i][j];
        end
    end
    ChebMatrix_inv = inv(ChebMatrix);
    return ChebMatrix_inv
end

function convert_matrix_to_interp(u,ppi::T where T<:Real,nodes_input)
    m = size(u,2);  one! = one(eltype(u))*one(eltype(nodes_input));
    k = size(u,1)-1;    ki = one!*k;
    ChebIdentity = Sequence(Chebyshev(k)^(k+1),one!*reshape(I[1:k+1,1:k+1],(k+1)^2));
    ChebMatrix_VectorOfVectors = coefficients.(ChebIdentity.(nodes_input));
    ChebMatrix = zeros(eltype(u),k+1,k+1)*one!;
    for j ∈ 1:k+1
        for i ∈ 1:k+1
            ChebMatrix[i,j] = ChebMatrix_VectorOfVectors[i][j];
        end
    end
    ChebMatrix_inv = inv(ChebMatrix);
    U_jℓ_entries = reshape(ChebMatrix_inv*u,(k+1)*m)
    U_jℓ = Sequence(Chebyshev(k)^m,U_jℓ_entries);
    return U_jℓ
end

function convert_matrix_to_interp(u,Cheb_matrix)
    m = size(u,2);
    k = size(u,1)-1;
    U_jℓ_entries = reshape(Cheb_matrix*u,(k+1)*m)
    U_jℓ = Sequence(Chebyshev(k)^m,U_jℓ_entries)
    return U_jℓ
end

function convert_DiffEqSol(a,k,m,δ;start=0,T=Float64::Type)
    mi = @interval m;   ppi = @interval pi;
    s = (0:m)/mi;       Δ = 1/mi;           ki = @interval k;
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    u_jℓ = zeros(T,k+2,m);           # Stored as prescribed type T.
    s_jℓ = zeros(Interval{T},k+2,m); # Stored as Interval of type T.
    for j ∈ 1:m
        for ℓ ∈ 1:k+2
            s_jℓ[ℓ,j] = s[j] + Δ*(nodes[ℓ]+1)/2
            u_jℓ[ℓ,j] = mid.(a(start + s_jℓ[ℓ,j]*δ)[1]);
        end
    end
    return s_jℓ,u_jℓ,reverse(nodes)
end

function make_s_grid(k,m,ppi;T=Float64::Type)
    mi = @interval m;   ki = @interval k;
    s = (0:m)/mi;       Δ = 1/mi;           
    s_jℓ = zeros(Interval{T},k+2,m); 
    nodes = reverse(cheb_nodes(k+1,ki+1,ppi));
    for j ∈ 1:m
        for ℓ ∈ 1:k+2
            s_jℓ[ℓ,j] = s[j] + Δ*(nodes[ℓ]+1)/2
        end
    end
    return s_jℓ
end

function num_components(c::Sequence{CartesianPower{Chebyshev}})
    base = length(component(c,1));
    len_c = length(coefficients(c))/base;
    return Int(len_c)
end

function dim_base(c::Sequence{CartesianPower{Chebyshev}})
    base = length(component(c,1));
    return Int(base)
end
