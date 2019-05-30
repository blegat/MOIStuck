module Bridges

include("../MathOptInterface.jl")
const MOI = MathOptInterface
const MOIU = MOI

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const SQF{T} = MOI.ScalarQuadraticFunction{T}
const VQF{T} = MOI.VectorQuadraticFunction{T}

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

abstract type AbstractBridge end
include("bridgeoptimizer.jl")
include("lazybridgeoptimizer.jl")

include("slackbridge.jl")
include("functionize_bridge.jl")
include("rsocbridge.jl")
include("socr_bridge.jl")
include("geomeanbridge.jl")
include("detbridge.jl")
include("soctopsdbridge.jl")

end # module
