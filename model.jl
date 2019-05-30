abstract type AbstractModel{T} <: MOI.ModelLike end
function MOI.add_constraint(model::AbstractModel, f::F, s::S) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    if MOI.supports_constraint(model, F, S)
        ci = MOI.ConstraintIndex{F, S}(model.nextconstraintid += 1)
        push!(model.constrmap, _add_constraint(model, ci, copy(f), copy(s)))
        return ci
    else
        throw(MOI.UnsupportedConstraint{F, S}())
    end
end

abstract type Constraints{F} end

struct ModelScalarConstraints{T, F <: (MOI).AbstractScalarFunction} <: Constraints{F}
end
function ModelScalarConstraints{T, F}() where {T, F}
    ModelScalarConstraints{T, F}()
end
const ConstraintEntry{F, S} = Tuple{MOI.ConstraintIndex{F, S}, F, S}
struct ModelVectorConstraints{T, F <: (MOI).AbstractVectorFunction} <: Constraints{F}
    moi_secondordercone::Vector{ConstraintEntry{F, MOI.SecondOrderCone}}
    moi_positivesemidefiniteconetriangle::Vector{ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}}
end
function ModelVectorConstraints{T, F}() where {T, F}
    ModelVectorConstraints{T, F}(ConstraintEntry{F, MOI.SecondOrderCone}[], ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}[])
end
mutable struct Model{T} <: AbstractModel{T}
    moi_vectoraffinefunction::ModelVectorConstraints{T, MOI.VectorAffineFunction{T}}
end
function Model{T}() where T
    Model{T}(
             ModelVectorConstraints{T, MOI.VectorAffineFunction{T}}()
            )
end
MOI.supports_constraint(model::Model{T}, ::Type{<:Union{}}, ::Type{<:Union{}}) where T = true
MOI.supports_constraint(model::Model{T}, ::Type{<:Union{MOI.VectorAffineFunction{T}}}, ::Type{<:Union{MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle}}) where T = true

function MOI.add_constraint(model::Model, f::F, s::S) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    if MOI.supports_constraint(model, F, S)
        ci = MOI.ConstraintIndex{F, S}(model.nextconstraintid += 1)
        push!(model.constrmap, _add_constraint(model, ci, copy(f), copy(s)))
        return ci
    else
        throw(MOI.UnsupportedConstraint{F, S}())
    end
end

function _add_constraint(constrs::Vector{ConstraintEntry{F, S}}, ci::MOI.ConstraintIndex, f::F,
                         s::S) where {F, S}
    push!(constrs, (ci, f, s))
    length(constrs)
end
function _add_constraint(model::ModelVectorConstraints, ci::(MOI.ConstraintIndex){F, <:MOI.SecondOrderCone}, args...) where F
    _add_constraint((model).moi_secondordercone, ci, args...)
end
function _add_constraint(model::ModelVectorConstraints, ci::(MOI.ConstraintIndex){F, <:MOI.PositiveSemidefiniteConeTriangle}, args...) where F
    _add_constraint((model).moi_positivesemidefiniteconetriangle, ci, args...)
end
function _add_constraint(model::Model, ci::(MOI.ConstraintIndex){<:MOI.VectorAffineFunction}, args...)
    _add_constraint((model).moi_vectoraffinefunction, ci, args...)
end
