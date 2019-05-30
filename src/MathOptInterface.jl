module MathOptInterface

abstract type ModelLike end
Base.broadcastable(model::ModelLike) = Ref(model)
abstract type AbstractOptimizer <: ModelLike end

abstract type UnsupportedError <: Exception end
abstract type NotAllowedError <: Exception end

struct ConstraintIndex{F, S}
    value::Int64
end
struct VariableIndex
    value::Int64
end

include("functions.jl")
include("sets.jl")

supports_constraint(model::ModelLike, ::Type{<:AbstractFunction}, ::Type{<:AbstractSet}) = false
struct UnsupportedConstraint{F<:AbstractFunction, S<:AbstractSet} <: UnsupportedError
    message::String # Human-friendly explanation why the attribute cannot be set
end
function add_constraint(model::ModelLike, func::AbstractFunction,
                        set::AbstractSet)
end

include("operate.jl")

end
