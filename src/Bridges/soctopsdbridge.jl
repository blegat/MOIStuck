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
_SOCtoPSDaff(f::MOI.VectorAffineFunction, ::Type) = _SOCtoPSDaff(f, MOIU.eachscalar(f)[1])

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
    f_scalars = MOIU.eachscalar(f)
    g = MOIU.operate!(*, T, f_scalars[2], convert(T, 2))
    _SOCtoPSDaff(f_scalars[[1; 3:n]], g)
end
