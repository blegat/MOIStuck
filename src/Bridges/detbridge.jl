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
