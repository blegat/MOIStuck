abstract type AbstractBridgeOptimizer <: MOI.AbstractOptimizer end

function is_bridged end
function is_bridged(b::AbstractBridgeOptimizer, ::Type{CI{F, S}}) where {F, S}
    return is_bridged(b, F, S)
end
is_bridged(b::AbstractBridgeOptimizer, ::Type{VI}) = false

function supports_bridging_constraint(::AbstractBridgeOptimizer,
                                      ::Type{<:MOI.AbstractFunction},
                                      ::Type{<:MOI.AbstractSet})
    return false
end

function bridge_type end

function concrete_bridge_type(b::AbstractBridgeOptimizer,
                              F::Type{<:MOI.AbstractFunction},
                              S::Type{<:MOI.AbstractSet})
    return concrete_bridge_type(bridge_type(b, F, S), F, S)
end

bridge(b::AbstractBridgeOptimizer, ci::CI) = b.bridges[ci.value]
function bridge(b::AbstractBridgeOptimizer,
                ci::CI{MOI.SingleVariable, S}) where S
    return b.single_variable_constraints[(ci.value, S)]
end

# Constraints
function MOI.supports_constraint(b::AbstractBridgeOptimizer,
                                F::Type{<:MOI.AbstractFunction},
                                S::Type{<:MOI.AbstractSet})
    if is_bridged(b, F, S)
        return supports_bridging_constraint(b, F, S)
    else
        return MOI.supports_constraint(b.model, F, S)
    end
end
function store_bridge(b::AbstractBridgeOptimizer, func::MOI.SingleVariable,
                      set::MOI.AbstractSet, bridge)
    b.single_variable_constraints[(func.variable.value, typeof(set))] = bridge
    return MOI.ConstraintIndex{MOI.SingleVariable, typeof(set)}(func.variable.value)
end
function store_bridge(b::AbstractBridgeOptimizer, func::MOI.AbstractFunction,
                      set::MOI.AbstractSet, bridge)
    push!(b.bridges, bridge)
    push!(b.constraint_types, (typeof(func), typeof(set)))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(b.bridges))
end
function MOI.add_constraint(b::AbstractBridgeOptimizer, f::MOI.AbstractFunction,
                            s::MOI.AbstractSet)
    println("Called")
    if is_bridged(b, typeof(f), typeof(s))
        # We compute `BridgeType` first as `concrete_bridge_type` calls
        # `bridge_type` which might throw an `UnsupportedConstraint` error in
        # which case, we do not want any modification to have been done
        BridgeType = concrete_bridge_type(b, typeof(f), typeof(s))
        # `add_constraint` might throw an `UnsupportedConstraint` but no
        # modification has been done in the previous line
        return store_bridge(b, f, s, bridge_constraint(BridgeType, b, f, s))
    else
        return MOI.add_constraint(b.model, f, s)
    end
end
function MOI.add_constraints(b::AbstractBridgeOptimizer, f::Vector{F},
                             s::Vector{S}) where { F <: MOI.AbstractFunction,
                             S <: MOI.AbstractSet}
    if is_bridged(b, F, S)
        return MOI.add_constraint.(b, f, s)
    else
        return MOI.add_constraints(b.model, f, s)
    end
end

# Variables
MOI.add_variable(b::AbstractBridgeOptimizer) = MOI.add_variable(b.model)
MOI.add_variables(b::AbstractBridgeOptimizer, n) = MOI.add_variables(b.model, n)
