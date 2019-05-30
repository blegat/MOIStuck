module Bridges

include("../MathOptInterface.jl")
const MOI = MathOptInterface

abstract type AbstractBridge end
abstract type AbstractBridgeOptimizer <: MOI.AbstractOptimizer end
function MOI.add_constraint(b::AbstractBridgeOptimizer, f::MOI.AbstractFunction,
                            s::MOI.AbstractSet)
    println("Called")
    if is_bridged(b, typeof(f), typeof(s))
        BridgeType = concrete_bridge_type(b, typeof(f), typeof(s))
        return store_bridge(b, f, s, bridge_constraint(BridgeType, b, f, s))
    else
        return MOI.add_constraint(b.model, f, s)
    end
end

mutable struct BridgeOptimizer{OT<:MOI.ModelLike} <: AbstractBridgeOptimizer
    model::OT   # Internal model
end
BridgeOptimizer(model::MOI.ModelLike) = BridgeOptimizer{typeof(model)}(model)
function is_bridged(b::BridgeOptimizer, F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    !MOI.supports_constraint(b.model, F, S)
end
function is_bridged(b::BridgeOptimizer, ::Type{MOI.ConstraintIndex{F, S}}) where {F, S}
    return is_bridged(b, F, S)
end
function concrete_bridge_type(b::AbstractBridgeOptimizer,
                              F::Type{<:MOI.AbstractFunction},
                              S::Type{<:MOI.AbstractSet})
    return concrete_bridge_type(bridge_type(b, F, S), F, S)
end
function MOI.supports_constraint(b::BridgeOptimizer,
                                F::Type{<:MOI.AbstractFunction},
                                S::Type{<:MOI.AbstractSet})
    if is_bridged(b, F, S)
        return supports_bridging_constraint(b, F, S)
    else
        return MOI.supports_constraint(b.model, F, S)
    end
end
function store_bridge(b::BridgeOptimizer, func::MOI.AbstractFunction,
                      set::MOI.AbstractSet, bridge)
    push!(b.bridges, bridge)
    push!(b.constraint_types, (typeof(func), typeof(set)))
    return MOI.ConstraintIndex{typeof(func), typeof(set)}(length(b.bridges))
end

const SVF = MOI.SingleVariable
const VVF = MOI.VectorOfVariables
const SAF{T} = MOI.ScalarAffineFunction{T}
const VAF{T} = MOI.VectorAffineFunction{T}
const SQF{T} = MOI.ScalarQuadraticFunction{T}
const VQF{T} = MOI.VectorQuadraticFunction{T}

const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

struct VectorSlackBridge{T, F, S} <: AbstractBridge
    slacks::Vector{MOI.VariableIndex}
    slacks_in_set::CI{MOI.VectorOfVariables, S}
    equality::CI{F, MOI.Zeros}
end
function bridge_constraint(::Type{VectorSlackBridge{T, F, S}}, model,
                           f::MOI.AbstractVectorFunction, s::S) where {T, F, S}
    d = MOI.dimension(s)
    slacks = MOI.VariableIndex.(1:d)
    new_f = MOI.operate(-, T, f, MOI.VectorOfVariables(slacks))
    slacks_in_set = MOI.add_constraint(model, MOI.VectorOfVariables(slacks), s)
    equality = MOI.add_constraint(model, new_f, MOI.Zeros(d))
    return VectorSlackBridge{T, F, S}(slacks, slacks_in_set, equality)
end
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = true
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.VectorOfVariables},
                        ::Type{<:MOI.Zeros}) where {T} = false
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.AbstractVectorFunction},
                        ::Type{<:MOI.Zeros}) where {T} = false
MOI.supports_constraint(::Type{VectorSlackBridge{T}},
                        ::Type{<:MOI.VectorOfVariables},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = false
function added_constraint_types(::Type{VectorSlackBridge{T, F, S}}) where {T, F<:MOI.AbstractVectorFunction, S}
    return [(F, MOI.Zeros), (MOI.VectorOfVariables, S)]
end
function concrete_bridge_type(::Type{<:VectorSlackBridge{T}},
                              F::Type{<:MOI.AbstractVectorFunction},
                              S::Type{<:MOI.AbstractVectorSet}) where T
    F2 = MOI.promote_operation(-, T, F, MOI.VectorOfVariables)
    return VectorSlackBridge{T, F2, S}
end


struct VectorFunctionizeBridge{T, S} <: AbstractBridge
    constraint::CI{MOI.VectorAffineFunction{T}, S}
end
function bridge_constraint(::Type{VectorFunctionizeBridge{T, S}}, model,
                           f::MOI.VectorOfVariables, s::S) where {T, S}
    constraint = MOI.add_constraint(model, MOI.VectorAffineFunction{T}(f), s)
    return VectorFunctionizeBridge{T, S}(constraint)
end

MOI.supports_constraint(::Type{VectorFunctionizeBridge{T}},
                        ::Type{MOI.VectorOfVariables},
                        ::Type{<:MOI.AbstractVectorSet}) where {T} = true
function added_constraint_types(::Type{VectorFunctionizeBridge{T, S}}) where {T, S}
    return [(MOI.VectorAffineFunction{T}, S)]
end
function concrete_bridge_type(::Type{<:VectorFunctionizeBridge{T}},
                              ::Type{<:MOI.AbstractVectorFunction},
                              S::Type{<:MOI.AbstractVectorSet}) where T
    return VectorFunctionizeBridge{T, S}
end


struct RSOCBridge{T, F, G} <: AbstractBridge
    soc::CI{F, MOI.SecondOrderCone}
end
function rotate_function(f::MOI.AbstractVectorFunction, T::Type)
    d = MOI.output_dimension(f)
    f_scalars = MOI.eachscalar(f)
    t = f_scalars[1]
    u = f_scalars[2]
    x = f_scalars[3:d]
    s2 = √2
    ts = MOI.operate!(/, T, t, s2)
    us = MOI.operate!(/, T, u, s2)
    # Cannot use `operate!` here since `ts` and `us` are needed for the next
    # line
    y  = ts - us
    z  = MOI.operate!(+, T, ts, us)
    return MOI.operate(vcat, T, z, y, x)
end
function bridge_constraint(::Type{RSOCBridge{T, F, G}}, model,
                           f::MOI.AbstractVectorFunction,
                           s::MOI.RotatedSecondOrderCone) where {T, F, G}
    soc = MOI.add_constraint(model, rotate_function(f, T),
                             MOI.SecondOrderCone(MOI.dimension(s)))
    return RSOCBridge{T, F, G}(soc)
end

function MOI.supports_constraint(::Type{RSOCBridge{T}},
                                ::Type{<:MOI.AbstractVectorFunction},
                                ::Type{MOI.RotatedSecondOrderCone}) where T
    return true
end
function added_constraint_types(::Type{<:RSOCBridge{T, F}}) where {T, F}
    return [(F, MOI.SecondOrderCone)]
end
function concrete_bridge_type(::Type{<:RSOCBridge{T}},
                              G::Type{<:MOI.AbstractVectorFunction},
                              ::Type{MOI.RotatedSecondOrderCone}) where T
    S = MOI.promote_operation(/, T, MOI.scalar_type(G), T)
    Y = MOI.promote_operation(-, T, S, S)
    Z = MOI.promote_operation(+, T, S, S)
    F = MOI.promote_operation(vcat, T, Z, Y, G)
    return RSOCBridge{T, F, G}
end


struct SOCRBridge{T, F, G} <: AbstractBridge
    rsoc::CI{F, MOI.RotatedSecondOrderCone}
end
function bridge_constraint(::Type{SOCRBridge{T, F, G}}, model,
                           f::MOI.AbstractVectorFunction,
                           s::MOI.SecondOrderCone) where {T, F, G}
    soc = MOI.add_constraint(model, rotate_function(f, T),
                             MOI.RotatedSecondOrderCone(MOI.dimension(s)))
    return SOCRBridge{T, F, G}(soc)
end

function MOI.supports_constraint(::Type{SOCRBridge{T}},
                                ::Type{<:MOI.AbstractVectorFunction},
                                ::Type{MOI.SecondOrderCone}) where T
    return true
end
function added_constraint_types(::Type{<:SOCRBridge{T, F}}) where {T, F}
    return [(F, MOI.RotatedSecondOrderCone)]
end
function concrete_bridge_type(::Type{<:SOCRBridge{T}},
                              G::Type{<:MOI.AbstractVectorFunction},
                              ::Type{MOI.SecondOrderCone}) where T
    S = MOI.promote_operation(/, T, MOI.scalar_type(G), T)
    Y = MOI.promote_operation(-, T, S, S)
    Z = MOI.promote_operation(+, T, S, S)
    F = MOI.promote_operation(vcat, T, Z, Y, G)
    return SOCRBridge{T, F, G}
end


function _ilog2(n, i)
    if n <= (one(n) << i)
        i
    else
        _ilog2(n, i+1)
    end
end
function ilog2(n::Integer)
    @assert n > zero(n)
    _ilog2(n, zero(n))
end

struct GeoMeanBridge{T, F, G} <: AbstractBridge
    # Initially, (t, x) is of dimension d so x is dimension (d-1)
    # We create n new variables so that there is 2^l = d-1+n variables x_i
    # We then need to create 2^l-1 new variables (1+2+...+2^{l-1})
    d::Int
    xij::Vector{VI}
    tubc::CI{F, MOI.LessThan{T}}
    socrc::Vector{CI{G, MOI.RotatedSecondOrderCone}}
end
function bridge_constraint(::Type{GeoMeanBridge{T, F, G}}, model,
                           f::MOI.AbstractVectorFunction,
                           s::MOI.GeometricMeanCone) where {T, F, G}
    d = s.dimension
    n = d-1
    l = ilog2(n)
    N = 1 << l
    xij = MOI.VariableIndex.(1:(N-1))
    f_scalars = MOI.eachscalar(f)

    xl1 = MOI.SingleVariable(xij[1])
    sN = one(T) / √N
    function _getx(i)
        if i > n
            return sN * xl1
        else
            return f_scalars[1+i]
        end
    end

    t = f_scalars[1]
    # With sqrt(2)^l*t - xl1, we should scale both the ConstraintPrimal and ConstraintDual
    tubc = MOI.add_scalar_constraint(model,
                                      MOI.operate!(+, T, t, -sN * xl1),
                                      MOI.LessThan(zero(T)),
                                      allow_modify_function=true)

    socrc = Vector{CI{G, MOI.RotatedSecondOrderCone}}(undef, N-1)
    offset = offsetnext = 0
    for i in 1:l
        offsetnext = offset + i
        for j in 1:(1 << (i-1))
            if i == l
                a = _getx(2j-1)
                b = _getx(2j)
            else
                a = one(T) * MOI.SingleVariable(xij[offsetnext+2j-1])
                b = one(T) * MOI.SingleVariable(xij[offsetnext+2j])
            end
            c = MOI.SingleVariable(xij[offset+j])
            socrc[offset + j] = MOI.add_constraint(model,
                                                   MOI.operate(vcat, T, a, b, c),
                                                   MOI.RotatedSecondOrderCone(3))
        end
        offset = offsetnext
    end
    GeoMeanBridge(d, xij, tubc, socrc)
end

function MOI.supports_constraint(::Type{GeoMeanBridge{T}},
                                ::Type{<:MOI.AbstractVectorFunction},
                                ::Type{MOI.GeometricMeanCone}) where T
    return true
end
function added_constraint_types(::Type{GeoMeanBridge{T, F, G}}) where {T, F, G}
    return [(F, MOI.LessThan{T}), (G, MOI.RotatedSecondOrderCone)]
end
function concrete_bridge_type(::Type{<:GeoMeanBridge{T}},
                              H::Type{<:MOI.AbstractVectorFunction},
                              ::Type{MOI.GeometricMeanCone}) where T
    S = MOI.scalar_type(H)
    A = MOI.promote_operation(*, T, T, MOI.SingleVariable)
    F = MOI.promote_operation(+, T, S, A)
    G = MOI.promote_operation(vcat, T, A, A, MOI.SingleVariable)
    return GeoMeanBridge{T, F, G}
end


struct RootDetBridge{T} <: AbstractBridge
    Δ::Vector{VI}
    sdindex::CI{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
    gmindex::CI{MOI.VectorAffineFunction{T}, MOI.GeometricMeanCone}
end
function bridge_constraint(::Type{RootDetBridge{T}}, model,
                           f::MOI.VectorOfVariables,
                           s::MOI.RootDetConeTriangle) where T
    return bridge_constraint(RootDetBridge{T}, model,
                             MOI.VectorAffineFunction{T}(f), s)
end
function bridge_constraint(::Type{RootDetBridge{T}}, model,
                           f::MOI.VectorAffineFunction{T},
                           s::MOI.RootDetConeTriangle) where T
    d = s.side_dimension
    tu, D, Δ, sdindex = extract_eigenvalues(model, f, d, 1)
    t = tu[1]
    DF = MOI.VectorAffineFunction{T}(MOI.VectorOfVariables(D))
    gmindex = MOI.add_constraint(model, MOI.operate(vcat, T, t, DF),
                                 MOI.GeometricMeanCone(d+1))

    return RootDetBridge(Δ, sdindex, gmindex)
end

MOI.supports_constraint(::Type{RootDetBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RootDetConeTriangle}) where T = true
added_constraint_types(::Type{RootDetBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RootDetConeTriangle}) where T = [(MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle), (MOI.VectorAffineFunction{T}, MOI.GeometricMeanCone)]



struct SOCtoPSDBridge{T} <: AbstractBridge
    dim::Int
    cr::CI{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end
function bridge_constraint(::Type{SOCtoPSDBridge{T}}, model, f,
                           s::MOI.SecondOrderCone) where T
    d = MOI.dimension(s)
    cr = MOI.add_constraint(model, _SOCtoPSDaff(f, T), MOI.PositiveSemidefiniteConeTriangle(d))
    SOCtoPSDBridge(d, cr)
end

_SOCtoPSDaff(f::MOI.VectorOfVariables, ::Type{T}) where T = _SOCtoPSDaff(MOI.VectorAffineFunction{T}(f), T)
_SOCtoPSDaff(f::MOI.VectorAffineFunction, ::Type) = _SOCtoPSDaff(f, MOI.eachscalar(f)[1])

MOI.supports_constraint(::Type{SOCtoPSDBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.SecondOrderCone}) where T = true
added_constraint_types(::Type{SOCtoPSDBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.SecondOrderCone}) where T = [(MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle)]

struct RSOCtoPSDBridge{T} <: AbstractBridge
    dim::Int
    cr::CI{MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle}
end

MOI.supports_constraint(::Type{RSOCtoPSDBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RotatedSecondOrderCone}) where T = true
added_constraint_types(::Type{RSOCtoPSDBridge{T}}, ::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}}, ::Type{MOI.RotatedSecondOrderCone}) where T = [(MOI.VectorAffineFunction{T}, MOI.PositiveSemidefiniteConeTriangle)]

function bridge_constraint(::Type{RSOCtoPSDBridge{T}}, model, f,
                           s::MOI.RotatedSecondOrderCone) where T
    d = MOI.dimension(s)-1
    cr = MOI.add_constraint(model, _RSOCtoPSDaff(f, T), MOI.PositiveSemidefiniteConeTriangle(d))
    RSOCtoPSDBridge(d, cr)
end

_RSOCtoPSDaff(f::MOI.VectorOfVariables, ::Type{T}) where T = _RSOCtoPSDaff(MOI.VectorAffineFunction{T}(f), T)
function _RSOCtoPSDaff(f::MOI.VectorAffineFunction, ::Type{T}) where T
    n = MOI.output_dimension(f)
    f_scalars = MOI.eachscalar(f)
    g = MOI.operate!(*, T, f_scalars[2], convert(T, 2))
    _SOCtoPSDaff(f_scalars[[1; 3:n]], g)
end

end # module
