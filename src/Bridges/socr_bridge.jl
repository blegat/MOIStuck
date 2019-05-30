struct SOCRBridge{T, F, G} <: AbstractBridge
    rsoc::CI{F, MOI.RotatedSecondOrderCone}
end
function bridge_constraint(::Type{SOCRBridge{T, F, G}}, model,
                           f::MOI.AbstractVectorFunction,
                           s::MOI.SecondOrderCone) where {T, F, G}
    soc = MOI.add_constraint(model, rotate_function(f, T),
                             MOI.RotatedSecondOrderCone(MOI.dimension(s)))
    return SOCRBridge{T, F, G}(soc)
end

function MOI.supports_constraint(::Type{SOCRBridge{T}},
                                ::Type{<:MOI.AbstractVectorFunction},
                                ::Type{MOI.SecondOrderCone}) where T
    return true
end
function added_constraint_types(::Type{<:SOCRBridge{T, F}}) where {T, F}
    return [(F, MOI.RotatedSecondOrderCone)]
end
function concrete_bridge_type(::Type{<:SOCRBridge{T}},
                              G::Type{<:MOI.AbstractVectorFunction},
                              ::Type{MOI.SecondOrderCone}) where T
    S = MOIU.promote_operation(/, T, MOIU.scalar_type(G), T)
    Y = MOIU.promote_operation(-, T, S, S)
    Z = MOIU.promote_operation(+, T, S, S)
    F = MOIU.promote_operation(vcat, T, Z, Y, G)
    return SOCRBridge{T, F, G}
end
