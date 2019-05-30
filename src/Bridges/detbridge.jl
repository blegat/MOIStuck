function trimap(i::Integer, j::Integer)
    if i < j
        trimap(j, i)
    else
        div((i-1)*i, 2) + j
    end
end

function extract_eigenvalues(model, f::MOI.VectorAffineFunction{T}, d::Int, offset::Int) where T
    f_scalars = MOIU.eachscalar(f)
    tu = [f_scalars[i] for i in 1:offset]

    n = trimap(d, d)
    X = f_scalars[offset .+ (1:n)]
    m = length(X.terms)
    M = m + n + d

    terms = Vector{MOI.VectorAffineTerm{T}}(undef, M)
    terms[1:m] = X.terms
    N = trimap(2d, 2d)
    constant = zeros(T, N); constant[1:n] = X.constants

    Δ = MOI.add_variables(model, n)

    cur = m
    for j in 1:d
        for i in j:d
            cur += 1
            terms[cur] = MOI.VectorAffineTerm(trimap(i, d + j),
                                              MOI.ScalarAffineTerm(one(T),
                                                                   Δ[trimap(i, j)]))
        end
        cur += 1
        terms[cur] = MOI.VectorAffineTerm(trimap(d + j, d + j),
                                          MOI.ScalarAffineTerm(one(T),
                                                               Δ[trimap(j, j)]))
    end
    @assert cur == M
    Y = MOI.VectorAffineFunction(terms, constant)
    sdindex = MOI.add_constraint(model, Y, MOI.PositiveSemidefiniteConeTriangle(2d))

    D = Δ[trimap.(1:d, 1:d)]

    return tu, D, Δ, sdindex
end

struct RootDetBridge{T} <: AbstractBridge
    Δ::Vector{VI}
    sdindex::CI{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
    gmindex::CI{MOI.VectorAffineFunction{T}, MOI.GeometricMeanCone}
end
function bridge_constraint(::Type{RootDetBridge{T}}, model,
                           f::MOI.VectorOfVariables,
                           s::MOI.RootDetConeTriangle) where T
    return bridge_constraint(RootDetBridge{T}, model,
                             MOI.VectorAffineFunction{T}(f), s)
end
function bridge_constraint(::Type{RootDetBridge{T}}, model,
                           f::MOI.VectorAffineFunction{T},
                           s::MOI.RootDetConeTriangle) where T
    d = s.side_dimension
    tu, D, Δ, sdindex = extract_eigenvalues(model, f, d, 1)
    t = tu[1]
    DF = MOI.VectorAffineFunction{T}(MOI.VectorOfVariables(D))
    gmindex = MOI.add_constraint(model, MOIU.operate(vcat, T, t, DF),
                                 MOI.GeometricMeanCone(d+1))

    return RootDetBridge(Δ, sdindex, gmindex)
end

MOI.supports_constraint(::Type{RootDetBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RootDetConeTriangle}) where T = true
added_constraint_types(::Type{RootDetBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RootDetConeTriangle}) where T = [(MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle), (MOI.VectorAffineFunction{T}, MOI.GeometricMeanCone)]
