struct VectorSlackBridge{T, F, S} <: AbstractBridge
    slacks::Vector{MOI.VariableIndex}
    slacks_in_set::CI{MOI.VectorOfVariables, S}
    equality::CI{F, MOI.Zeros}
end
function bridge_constraint(::Type{VectorSlackBridge{T, F, S}}, model,
                           f::MOI.AbstractVectorFunction, s::S) where {T, F, S}
    d = MOI.dimension(s)
    slacks = MOI.add_variables(model, d)
    new_f = MOIU.operate(-, T, f, MOI.VectorOfVariables(slacks))
    slacks_in_set = MOI.add_constraint(model, MOI.VectorOfVariables(slacks), s)
    equality = MOI.add_constraint(model, new_f, MOI.Zeros(d))
    return VectorSlackBridge{T, F, S}(slacks, slacks_in_set, equality)
end

MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = true
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.VectorOfVariables},
                        ::Type{<:MOI.Zeros}) where {T} = false
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:MOI.Zeros}) where {T} = false
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.VectorOfVariables},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = false
function added_constraint_types(::Type{VectorSlackBridge{T, F, S}}) where {T, F<:MOI.AbstractVectorFunction, S}
    return [(F, MOI.Zeros), (MOI.VectorOfVariables, S)]
end
function concrete_bridge_type(::Type{<:VectorSlackBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              S::Type{<:MOI.AbstractVectorSet}) where T
    F2 = MOIU.promote_operation(-, T, F, MOI.VectorOfVariables)
    return VectorSlackBridge{T, F2, S}
end
