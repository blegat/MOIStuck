struct VectorFunctionizeBridge{T, S} <: AbstractBridge
    constraint::CI{MOI.VectorAffineFunction{T}, S}
end
function bridge_constraint(::Type{VectorFunctionizeBridge{T, S}}, model,
                           f::MOI.VectorOfVariables, s::S) where {T, S}
    constraint = MOI.add_constraint(model, MOI.VectorAffineFunction{T}(f), s)
    return VectorFunctionizeBridge{T, S}(constraint)
end

MOI.supports_constraint(::Type{VectorFunctionizeBridge{T}},
                        ::Type{MOI.VectorOfVariables},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = true
function added_constraint_types(::Type{VectorFunctionizeBridge{T, S}}) where {T, S}
    return [(MOI.VectorAffineFunction{T}, S)]
end
function concrete_bridge_type(::Type{<:VectorFunctionizeBridge{T}},
                              ::Type{<:MOI.AbstractVectorFunction},
                              S::Type{<:MOI.AbstractVectorSet}) where T
    return VectorFunctionizeBridge{T, S}
end
