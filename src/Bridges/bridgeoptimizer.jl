"""
    AbstractBridgeOptimizer

A bridge optimizer applies given constraint bridges to a given optimizer thus
extending the types of supported constraints. The attributes of the inner
optimizer are automatically transformed to make the bridges transparent, e.g.
the variables and constraints created by the bridges are hidden.

By convention, the inner optimizer should be stored in a `model` field and
the dictionary mapping constraint indices to bridges should be stored in a
`bridges` field. If a bridge optimizer deviates from these conventions, it
should implement the functions `MOI.optimize!` and `bridge` respectively.
"""
abstract type AbstractBridgeOptimizer <: MOI.AbstractOptimizer end

# AbstractBridgeOptimizer interface

"""
    is_bridged(b::AbstractBridgeOptimizer, F::Type{<:MOI.AbstractFunction},
              S::Type{<:MOI.AbstractSet})::Bool

Return a `Bool` indicating whether `b` tries to bridge `F`-in-`S` constraints
instead of passing it as is to its internal model.
"""
function is_bridged end
# Syntactic sugar
function is_bridged(b::AbstractBridgeOptimizer, ::Type{CI{F, S}}) where {F, S}
    return is_bridged(b, F, S)
end
# We don't bridge variables.
is_bridged(b::AbstractBridgeOptimizer, ::Type{VI}) = false

"""
    supports_bridging_constraint(b::AbstractBridgeOptimizer,
                               F::Type{<:MOI.AbstractFunction},
                               S::Type{<:MOI.AbstractSet})::Bool

Return a `Bool` indicating whether `b` supports bridging `F`-in-`S` constraints.
"""
function supports_bridging_constraint(::AbstractBridgeOptimizer,
                                      ::Type{<:MOI.AbstractFunction},
                                      ::Type{<:MOI.AbstractSet})
    return false
end

"""
    bridge_type(b::AbstractBridgeOptimizer,
                F::Type{<:MOI.AbstractFunction},
                S::Type{<:MOI.AbstractSet})

Return the `AbstractBridge` type to be used to bridge `F`-in-`S` constraints.
This function should only be called if `is_bridged(b, F, S)`.
"""
function bridge_type end

function concrete_bridge_type(b::AbstractBridgeOptimizer,
                              F::Type{<:MOI.AbstractFunction},
                              S::Type{<:MOI.AbstractSet})
    return concrete_bridge_type(bridge_type(b, F, S), F, S)
end

"""
    bridge(b::AbstractBridgeOptimizer, ci::CI)

Return the `AbstractBridge` used to bridge the constraint with index `ci`.
"""
bridge(b::AbstractBridgeOptimizer, ci::CI) = b.bridges[ci.value]
# By convention, they should be stored in a `bridges` field using a
# dictionary-like object.
function bridge(b::AbstractBridgeOptimizer,
                ci::CI{MOI.SingleVariable, S}) where S
    return b.single_variable_constraints[(ci.value, S)]
end

# Implementation of the MOI interface for AbstractBridgeOptimizer

MOI.optimize!(b::AbstractBridgeOptimizer) = MOI.optimize!(b.model)
# By convention, the model should be stored in a `model` field

function MOI.is_empty(b::AbstractBridgeOptimizer)
    return isempty(b.bridges) && isempty(b.constraint_types) &&
           MOI.is_empty(b.model) && isempty(b.single_variable_constraints)
end
function MOI.empty!(b::AbstractBridgeOptimizer)
    MOI.empty!(b.model)
    empty!(b.bridges)
    empty!(b.constraint_types)
    empty!(b.single_variable_constraints)
    empty!(b.con_to_name)
    b.name_to_con = nothing
end
function MOI.supports(b::AbstractBridgeOptimizer,
                      attr::Union{MOI.AbstractModelAttribute,
                                  MOI.AbstractOptimizerAttribute})
    return MOI.supports(b.model, attr)
end

function MOI.copy_to(mock::AbstractBridgeOptimizer, src::MOI.ModelLike; kws...)
    MOIU.automatic_copy_to(mock, src; kws...)
end

# References
function _constraint_index(b::AbstractBridgeOptimizer, i::Integer)
    F, S = b.constraint_types[i]
    return MOI.ConstraintIndex{F, S}(i)
end
MOI.is_valid(b::AbstractBridgeOptimizer, vi::VI) = MOI.is_valid(b.model, vi)
function MOI.is_valid(b::AbstractBridgeOptimizer, ci::CI{F, S}) where {F, S}
    if is_bridged(b, typeof(ci))
        return 1 ≤ ci.value ≤ length(b.bridges) && b.bridges[ci.value] !== nothing &&
            (F, S) == b.constraint_types[ci.value]
    else
        return MOI.is_valid(b.model, ci)
    end
end
function MOI.is_valid(b::AbstractBridgeOptimizer,
                      ci::CI{MOI.SingleVariable, S}) where S
    if is_bridged(b, typeof(ci))
        return haskey(b.single_variable_constraints, (ci.value, S))
    else
        return MOI.is_valid(b.model, ci)
    end
end
function MOI.delete(b::AbstractBridgeOptimizer, vi::VI)
    # Delete all `MOI.SingleVariable` constraint of this variable
    for key in keys(b.single_variable_constraints)
        if key[1] == vi.value
            MOI.delete(b, MOI.ConstraintIndex{MOI.SingleVariable,
                                              key[2]}(vi.value))
        end
    end
    MOI.delete(b.model, vi)
end
function _remove_bridge(b, ci)
    b.bridges[ci.value] = nothing
end
function _remove_bridge(b, ci::MOI.ConstraintIndex{MOI.SingleVariable, S}) where S
    delete!(b.single_variable_constraints, (ci.value, S))
end
function MOI.delete(b::AbstractBridgeOptimizer, ci::CI)
    if is_bridged(b, typeof(ci))
        if !MOI.is_valid(b, ci)
            throw(MOI.InvalidIndex(ci))
        end
        MOI.delete(b, bridge(b, ci))
        _remove_bridge(b, ci)
        b.name_to_con = nothing
        if haskey(b.con_to_name, ci)
            delete!(b.con_to_name, ci)
        end
    else
        MOI.delete(b.model, ci)
    end
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
function MOI.modify(b::AbstractBridgeOptimizer, ci::CI,
                     change::MOI.AbstractFunctionModification)
    if is_bridged(b, typeof(ci))
        MOI.modify(b, bridge(b, ci), change)
    else
        MOI.modify(b.model, ci, change)
    end
end

# Objective
function MOI.modify(b::AbstractBridgeOptimizer, obj::MOI.ObjectiveFunction,
                     change::MOI.AbstractFunctionModification)
    MOI.modify(b.model, obj, change)
end

# Variables
MOI.add_variable(b::AbstractBridgeOptimizer) = MOI.add_variable(b.model)
MOI.add_variables(b::AbstractBridgeOptimizer, n) = MOI.add_variables(b.model, n)

# TODO add transform
