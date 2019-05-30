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
    f_scalars = MOIU.eachscalar(f)

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
    tubc = MOIU.add_scalar_constraint(model,
                                      MOIU.operate!(+, T, t, -sN * xl1),
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
                                                   MOIU.operate(vcat, T, a, b, c),
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
    S = MOIU.scalar_type(H)
    A = MOIU.promote_operation(*, T, T, MOI.SingleVariable)
    F = MOIU.promote_operation(+, T, S, A)
    G = MOIU.promote_operation(vcat, T, A, A, MOI.SingleVariable)
    return GeoMeanBridge{T, F, G}
end
