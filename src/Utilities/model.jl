const ConstraintEntry{F, S} = Tuple{CI{F, S}, F, S}

function _add_constraint(constrs::Vector{ConstraintEntry{F, S}}, ci::CI, f::F,
                         s::S) where {F, S}
    push!(constrs, (ci, f, s))
    length(constrs)
end

abstract type AbstractModel{T} <: MOI.ModelLike end

function MOI.add_constraint(model::AbstractModel, f::F, s::S) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    if MOI.supports_constraint(model, F, S)
        ci = CI{F, S}(model.nextconstraintid += 1)
        push!(model.constrmap, _add_constraint(model, ci, copy(f), copy(s)))
        return ci
    else
        throw(MOI.UnsupportedConstraint{F, S}())
    end
end

abstract type Constraints{F} end
