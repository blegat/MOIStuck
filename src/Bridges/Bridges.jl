module Bridges

include("../Utilities/Utilities.jl")
const MOIU = Utilities
const MOI = Utilities.MathOptInterface

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const SQF{T} = MOI.ScalarQuadraticFunction{T}
const VQF{T} = MOI.VectorQuadraticFunction{T}

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

include("bridge.jl")
include("bridgeoptimizer.jl")
include("lazybridgeoptimizer.jl")

include("flip_sign_bridge.jl")
include("vectorizebridge.jl")
include("scalarizebridge.jl")
include("slackbridge.jl")
include("functionize_bridge.jl")
include("intervalbridge.jl")
include("rsocbridge.jl")
include("socr_bridge.jl")
include("quadtosocbridge.jl")
include("geomeanbridge.jl")
include("square_bridge.jl")
include("detbridge.jl")
include("soctopsdbridge.jl")

end # module
