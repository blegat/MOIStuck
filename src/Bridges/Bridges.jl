module Bridges

include("../MathOptInterface.jl")
const MOI = MathOptInterface
const MOIU = MOI

abstract type AbstractBridge end
abstract type AbstractBridgeOptimizer <: MOI.AbstractOptimizer end
function MOI.add_constraint(b::AbstractBridgeOptimizer, f::MOI.AbstractFunction,
                            s::MOI.AbstractSet)
    println("Called")
    if is_bridged(b, typeof(f), typeof(s))
        BridgeType = concrete_bridge_type(b, typeof(f), typeof(s))
        return store_bridge(b, f, s, bridge_constraint(BridgeType, b, f, s))
    else
        return MOI.add_constraint(b.model, f, s)
    end
end

mutable struct BridgeOptimizer{OT<:MOI.ModelLike} <: AbstractBridgeOptimizer
    model::OT   # Internal model
end
BridgeOptimizer(model::MOI.ModelLike) = BridgeOptimizer{typeof(model)}(model)
function is_bridged(b::BridgeOptimizer, F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    !MOI.supports_constraint(b.model, F, S)
end
function is_bridged(b::BridgeOptimizer, ::Type{MOI.ConstraintIndex{F, S}}) where {F, S}
    return is_bridged(b, F, S)
end
function concrete_bridge_type(b::AbstractBridgeOptimizer,
                              F::Type{<:MOI.AbstractFunction},
                              S::Type{<:MOI.AbstractSet})
    return concrete_bridge_type(bridge_type(b, F, S), F, S)
end
function MOI.supports_constraint(b::BridgeOptimizer,
                                F::Type{<:MOI.AbstractFunction},
                                S::Type{<:MOI.AbstractSet})
    if is_bridged(b, F, S)
        return supports_bridging_constraint(b, F, S)
    else
        return MOI.supports_constraint(b.model, F, S)
    end
end
function store_bridge(b::BridgeOptimizer, func::MOI.AbstractFunction,
                      set::MOI.AbstractSet, bridge)
    push!(b.bridges, bridge)
    push!(b.constraint_types, (typeof(func), typeof(set)))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(b.bridges))
end

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const SQF{T} = MOI.ScalarQuadraticFunction{T}
const VQF{T} = MOI.VectorQuadraticFunction{T}

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

include("slackbridge.jl")
include("functionize_bridge.jl")
include("rsocbridge.jl")
include("socr_bridge.jl")
include("geomeanbridge.jl")
include("detbridge.jl")
include("soctopsdbridge.jl")

end # module
