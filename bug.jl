include("src/Bridges/Bridges.jl")
const MOIB = Bridges
const MOIU = MOIB.Utilities
const MOI = MOIU.MathOptInterface
struct ModelScalarConstraints{T, F <: (MOI).AbstractScalarFunction} <: MOIU.Constraints{F}
end
function ModelScalarConstraints{T, F}() where {T, F}
    ModelScalarConstraints{T, F}()
end
struct ModelVectorConstraints{T, F <: (MOI).AbstractVectorFunction} <: MOIU.Constraints{F}
    moi_secondordercone::MOIU.Vector{MOIU.ConstraintEntry{F, MOI.SecondOrderCone}}
    moi_positivesemidefiniteconetriangle::MOIU.Vector{MOIU.ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}}
end
function ModelVectorConstraints{T, F}() where {T, F}
    ModelVectorConstraints{T, F}(MOIU.ConstraintEntry{F, MOI.SecondOrderCone}[], MOIU.ConstraintEntry{F, MOI.PositiveSemidefiniteConeTriangle}[])
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
