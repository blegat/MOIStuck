include("src/Bridges/Bridges.jl")
const MOIB = Bridges
const MOIU = MOIB.Utilities
const MOI = MOIU.MathOptInterface
MOIU.@model(NoRSOCModel, (), (),
            (MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle),
            (), (), (), (),
            (MOI.VectorAffineFunction,))
mock = NoRSOCModel{Float64}()
bridged_mock = MOIB.LazyBridgeOptimizer(mock)
tQ = MOI.add_variables(bridged_mock, 4)
vov = MOI.VectorOfVariables(tQ)
println("Calling")
cX = MOI.add_constraint(bridged_mock, vov, MOI.RootDetConeTriangle(2))
