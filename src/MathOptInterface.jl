module MathOptInterface

abstract type ModelLike end
Base.broadcastable(model::ModelLike) = Ref(model)
abstract type AbstractOptimizer <: ModelLike end

abstract type UnsupportedError <: Exception end
abstract type NotAllowedError <: Exception end

include("indextypes.jl")
include("functions.jl")
include("sets.jl")
include("constraints.jl")

include("Utilities/Utilities.jl")

end
