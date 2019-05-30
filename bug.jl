include("src/Bridges/Bridges.jl")
const MOIB = Bridges
const MOIU = MOIB.Utilities
const MOI = MOIU.MathOptInterface
struct ModelScalarConstraints{T, F <: (MOI).AbstractScalarFunction} <: MOIU.Constraints{F}
end
function ModelScalarConstraints{T, F}() where {T, F}
    ModelScalarConstraints{T, F}()
end
const ConstraintEntry{F, S} = Tuple{MOI.ConstraintIndex{F, S}, F, S}
struct ModelVectorConstraints{T, F <: (MOI).AbstractVectorFunction} <: MOIU.Constraints{F}
    moi_secondordercone::MOIU.Vector{ConstraintEntry{F, MOI.SecondOrderCone}}
    moi_positivesemidefiniteconetriangle::MOIU.Vector{ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}}
end
function ModelVectorConstraints{T, F}() where {T, F}
    ModelVectorConstraints{T, F}(ConstraintEntry{F, MOI.SecondOrderCone}[], ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}[])
end
mutable struct Model{T} <: MOIU.AbstractModel{T}
    moi_vectoraffinefunction::ModelVectorConstraints{T, MOI.VectorAffineFunction{T}}
end
function Model{T}() where T
    Model{T}(
             ModelVectorConstraints{T, MOI.VectorAffineFunction{T}}()
            )
end
MOI.supports_constraint(model::Model{T}, ::MOIU.Type{<:MOIU.Union{}}, ::MOIU.Type{<:MOIU.Union{}}) where T = true
MOI.supports_constraint(model::Model{T}, ::MOIU.Type{<:MOIU.Union{MOI.VectorAffineFunction{T}}}, ::MOIU.Type{<:MOIU.Union{MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle}}) where T = true

function MOI.add_variable(model::Model{T}) where T
    vi = VI(model.num_variables_created += 1)
    push!(model.single_variable_mask, 0x0)
    push!(model.lower_bound, zero(T))
    push!(model.upper_bound, zero(T))
    if model.variable_indices !== nothing
        push!(model.variable_indices, vi)
    end
    return vi
end
function MOI.add_variables(model::Model, n::Integer)
    return [MOI.add_variable(model) for i in 1:n]
end

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
function MOIU._add_constraint(model::ModelVectorConstraints, ci::(MOI.ConstraintIndex){F, <:MOI.SecondOrderCone}, args...) where F
    MOIU._add_constraint((model).moi_secondordercone, ci, args...)
end
function MOIU._add_constraint(model::ModelVectorConstraints, ci::(MOI.ConstraintIndex){F, <:MOI.PositiveSemidefiniteConeTriangle}, args...) where F
    MOIU._add_constraint((model).moi_positivesemidefiniteconetriangle, ci, args...)
end
function MOIU._add_constraint(model::Model, ci::(MOI.ConstraintIndex){<:MOI.VectorAffineFunction}, args...)
    MOIU._add_constraint((model).moi_vectoraffinefunction, ci, args...)
end
mock = Model{Float64}()
bridged_mock = MOIB.LazyBridgeOptimizer(mock)
vov = MOI.VectorOfVariables(MOI.VariableIndex.(1:4))
println("Calling")
cX = MOI.add_constraint(bridged_mock, vov, MOI.RootDetConeTriangle(2))
