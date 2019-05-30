include("src/Bridges/Bridges.jl")
const MOIB = Bridges
const MOIU = MOIB.Utilities
const MOI = MOIU.MathOptInterface
include("model.jl")
mock = Model{Float64}()
bridged_mock = MOIB.LazyBridgeOptimizer(mock)
vov = MOI.VectorOfVariables(MOI.VariableIndex.(1:4))
println("Calling")
cX = MOI.add_constraint(bridged_mock, vov, MOI.RootDetConeTriangle(2))
