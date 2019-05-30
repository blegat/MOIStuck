# Functions convertible to a ScalarAffineFunction
const ScalarAffineLike{T} = Union{T, MOI.SingleVariable, MOI.ScalarAffineFunction{T}}
# Functions convertible to a ScalarQuadraticFunction
const ScalarQuadraticLike{T} = Union{ScalarAffineLike{T}, MOI.ScalarQuadraticFunction{T}}

# Used for overloading Base operator functions so `T` is not in the union to
# avoid overloading e.g. `+(::Float64, ::Float64)`
const ScalarLike{T} = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{T},
                            MOI.ScalarQuadraticFunction{T}}

# Functions convertible to a VectorAffineFunction
const VectorAffineLike{T} = Union{Vector{T}, MOI.VectorOfVariables, MOI.VectorAffineFunction{T}}
# Functions convertible to a VectorQuadraticFunction
const VectorQuadraticLike{T} = Union{VectorAffineLike{T}, MOI.VectorQuadraticFunction{T}}

# Used for overloading Base operator functions so `T` is not in the union to
# avoid overloading e.g. `+(::Float64, ::Float64)`
const VectorLike{T} = Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{T},
                            MOI.VectorQuadraticFunction{T}}

###################################### +/- #####################################
## promote_operation

function promote_operation(::typeof(-), ::Type{T},
                           ::Type{<:ScalarAffineLike{T}}) where T
    return MOI.ScalarAffineFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:ScalarAffineLike{T}},
                           ::Type{<:ScalarAffineLike{T}}) where T
    return MOI.ScalarAffineFunction{T}
end
function promote_operation(::typeof(-), ::Type{T},
                           ::Type{<:ScalarQuadraticLike{T}}) where T
    return MOI.ScalarQuadraticFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:ScalarQuadraticLike{T}},
                           ::Type{<:ScalarQuadraticLike{T}}) where T
    return MOI.ScalarQuadraticFunction{T}
end

## operate!
# + with at least 3 arguments
function operate!(op::typeof(+), ::Type{T}, f, g, h, args...) where T
    operate!(op, T, f, g)
    return operate!(+, T, f, h, args...)
end

# Unary -
function operate!(op::typeof(-), ::Type{T}, f::MOI.ScalarQuadraticFunction{T}) where T
    operate_terms!(-, f.quadratic_terms)
    operate_terms!(-, f.affine_terms)
    f.constant = -f.constant
    return f
end


# Scalar Variable +/- ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.SingleVariable,
                  g::ScalarQuadraticLike) where T
    return operate(op, T, f, g)
end
# Scalar Affine +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarAffineFunction{T},
                  g::T) where T
    f.constant = op(f.constant, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarAffineFunction{T},
                  g::MOI.SingleVariable) where T
    push!(f.terms, MOI.ScalarAffineTerm(op(one(T)), g.variable))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarAffineFunction{T},
                  g::MOI.ScalarAffineFunction{T}) where T
    append!(f.terms, operate_terms(op, g.terms))
    f.constant = op(f.constant, g.constant)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarAffineFunction{T},
                  g::MOI.ScalarQuadraticFunction{T}) where T
    return operate(op, T, f, g)
end
# Scalar Quadratic +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarQuadraticFunction{T},
                  g::T) where T
    f.constant = op(f.constant, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarQuadraticFunction{T},
                  g::MOI.SingleVariable) where T
    push!(f.affine_terms, MOI.ScalarAffineTerm(op(one(T)), g.variable))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarQuadraticFunction{T},
                  g::MOI.ScalarAffineFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.terms))
    f.constant = op(f.constant, g.constant)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.ScalarQuadraticFunction{T},
                  g::MOI.ScalarQuadraticFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.affine_terms))
    append!(f.quadratic_terms, operate_terms(op, g.quadratic_terms))
    f.constant = op(f.constant, g.constant)
    return f
end

## operate
# + with at least 3 arguments, can use in-place as the user cannot use
# intermediate results
function operate(op::typeof(+), ::Type{T}, f, g, h, args...) where T
    return operate!(+, T, operate(+, T, f, g), h, args...)
end

# Unary +
function operate(op::typeof(+), ::Type{T}, f::MOI.AbstractFunction) where T
    return f
end

# Scalar number +/- ...
function operate(op::typeof(+), ::Type{T}, α::T, f::ScalarLike{T}) where T
    return operate(op, T, f, α)
end
function operate(op::typeof(-), ::Type{T}, α::T, f::ScalarLike{T}) where T
    return operate!(+, T, operate(-, T, f), α)
end

# Scalar Variable +/- ...
function operate(::typeof(-), ::Type{T}, f::MOI.SingleVariable) where T
    return MOI.ScalarAffineFunction{T}(
        [MOI.ScalarAffineTerm(-one(T), f.variable)], zero(T))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.SingleVariable, α::T) where T
    return MOI.ScalarAffineFunction{T}(
        [MOI.ScalarAffineTerm(one(T), f.variable)], op(α))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.SingleVariable,
                 g::MOI.SingleVariable) where T
    return MOI.ScalarAffineFunction{T}(
        [MOI.ScalarAffineTerm(one(T), f.variable),
         MOI.ScalarAffineTerm(op(one(T)), g.variable)],
        zero(T))
end
function operate(op::typeof(+), ::Type{T},
                 f::MOI.SingleVariable,
                 g::Union{MOI.ScalarAffineFunction{T},
                          MOI.ScalarQuadraticFunction{T}}) where T
    return operate(op, T, g, f)
end
function operate(op::typeof(-), ::Type{T},
                 f::MOI.SingleVariable,
                 g::Union{MOI.ScalarAffineFunction{T},
                          MOI.ScalarQuadraticFunction{T}}) where T
    return operate!(+, T, operate(-, T, g), f)
end
# Scalar Affine +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
                 f::MOI.ScalarAffineFunction{T}) where T
    return MOI.ScalarAffineFunction(operate_terms(op, f.terms),
                                    op(f.constant))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.ScalarAffineFunction{T},
                 g::ScalarAffineLike{T}) where T
    return operate!(op, T, copy(f), g)
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.ScalarAffineFunction{T},
                 g::MOI.ScalarQuadraticFunction{T}) where T
    MOI.ScalarQuadraticFunction([f.terms; operate_terms(op, g.affine_terms)],
                                operate_terms(op, g.quadratic_terms),
                                op(f.constant, g.constant))
end
# Scalar Quadratic +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
                 f::MOI.ScalarQuadraticFunction{T}) where T
    return MOI.ScalarQuadraticFunction(
        operate_terms(op, f.affine_terms),
        operate_terms(op, f.quadratic_terms),
        op(f.constant))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.ScalarQuadraticFunction{T},
                 g::ScalarQuadraticLike{T}) where T
    operate!(op, T, copy(f), g)
end

function Base.:+(args::ScalarLike{T}...) where T
    return operate(+, T, args...)
end
function Base.:+(α::T, f::ScalarLike{T}...) where T
    return operate(+, T, α, f...)
end
function Base.:+(f::ScalarLike{T}, α::T) where T
    return operate(+, T, f, α)
end
function Base.:-(args::ScalarLike{T}...) where T
    return operate(-, T, args...)
end
function Base.:-(f::ScalarLike{T}, α::T) where T
    return operate(-, T, f, α)
end
function Base.:-(α::T, f::ScalarLike{T}) where T
    return operate(-, T, α, f)
end

# Vector +/-
###############################################################################
function promote_operation(::typeof(-), ::Type{T},
    ::Type{<:VectorAffineLike{T}}) where T
    return MOI.VectorAffineFunction{T}
end
function promote_operation(::typeof(-), ::Type{T},
    ::Type{<:VectorQuadraticLike{T}}) where T
    return MOI.VectorQuadraticFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:VectorAffineLike{T}},
                           ::Type{<:VectorAffineLike{T}}) where T
    return MOI.VectorAffineFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:VectorQuadraticLike{T}},
                           ::Type{<:VectorQuadraticLike{T}}) where T
    return MOI.VectorQuadraticFunction{T}
end

# Vector Variable +/- ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorOfVariables,
                  g::VectorQuadraticLike) where T
    return operate(op, T, f, g)
end
# Vector Affine +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorAffineFunction{T},
                  g::Vector{T}) where T
    @assert MOI.output_dimension(f) == length(g)
    f.constants .= op.(f.constants, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorAffineFunction{T},
                  g::MOI.VectorOfVariables) where T
    d = MOI.output_dimension(g)
    @assert MOI.output_dimension(f) == d
    append!(f.terms, MOI.VectorAffineTerm.(
                         collect(1:d),
                         MOI.ScalarAffineTerm.(op(one(T)), g.variables)))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorAffineFunction{T},
                  g::MOI.VectorAffineFunction{T}) where T
    append!(f.terms, operate_terms(op, g.terms))
    f.constants .= op.(f.constants, g.constants)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorAffineFunction{T},
                  g::MOI.VectorQuadraticFunction{T}) where T
    return operate(op, T, f, g)
end
# Vector Quadratic +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorQuadraticFunction{T},
                  g::Vector{T}) where T
    @assert MOI.output_dimension(f) == length(g)
    f.constants .= op.(f.constants, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorQuadraticFunction{T},
                  g::MOI.VectorOfVariables) where T
    d = MOI.output_dimension(g)
    @assert MOI.output_dimension(f) == d
    append!(f.affine_terms, MOI.VectorAffineTerm.(collect(1:d), MOI.ScalarAffineTerm.(op(one(T)), g.variables)))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorQuadraticFunction{T},
                  g::MOI.VectorAffineFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.terms))
    f.constants .= op.(f.constants, g.constants)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::MOI.VectorQuadraticFunction{T},
                  g::MOI.VectorQuadraticFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.affine_terms))
    append!(f.quadratic_terms, operate_terms(op, g.quadratic_terms))
    f.constants .= op.(f.constants, g.constants)
    return f
end

## operate
# + with at least 3 arguments, can use in-place as the user cannot use
# intermediate results
# overload
# function operate(op::typeof(+), ::Type{T}, f, g, h, args...) where T
#     return operate!(+, T, operate(+, T, f, g), h, args...)
# end

# function operate(op::typeof(+), ::Type{T}, f::VectorOfVariables) where T
#     return f
# end

function operate(op::typeof(-), ::Type{T}, f::MOI.VectorOfVariables) where T
    d = MOI.output_dimension(f)
    return MOI.VectorAffineFunction{T}(
               MOI.VectorAffineTerm.(
                   collect(1:d),
                   MOI.ScalarAffineTerm.(-one(T), f.variables)),
               fill(zero(T),d))
end

# Vector number +/- ...
function operate(op::typeof(+), ::Type{T}, α::Vector{T}, f::VectorLike{T}) where T
    return operate(op, T, f, α)
end
function operate(op::typeof(-), ::Type{T}, α::Vector{T}, f::VectorLike{T}) where T
    return operate!(+, T, operate(-, T, f), α)
end

# Vector Variable +/- ...
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::MOI.VectorOfVariables, α::Vector{T}) where T
    d = MOI.output_dimension(f)
    @assert length(α) == d
    return MOI.VectorAffineFunction{T}(
               MOI.VectorAffineTerm.(
                   collect(1:d),
                   MOI.ScalarAffineTerm.(one(T), f.variables)),
               op.(α))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::MOI.VectorOfVariables,
                 g::MOI.VectorOfVariables) where T
    d = MOI.output_dimension(f)
    @assert MOI.output_dimension(g) == d
    return MOI.VectorAffineFunction{T}(
               vcat(
                   MOI.VectorAffineTerm.(
                       collect(1:d),
                       MOI.ScalarAffineTerm.(one(T), f.variables)),
                   MOI.VectorAffineTerm.(
                       collect(1:d),
                       MOI.ScalarAffineTerm.(op(one(T)), g.variables))),
               fill(zero(T),d))
end
function operate(op::typeof(+), ::Type{T},
    f::MOI.VectorOfVariables,
    g::Union{MOI.VectorAffineFunction{T},
                MOI.VectorQuadraticFunction{T}}) where T
    return operate(op, T, g, f)
end
function operate(op::typeof(-), ::Type{T},
    f::MOI.VectorOfVariables,
    g::Union{MOI.VectorAffineFunction{T},
                MOI.VectorQuadraticFunction{T}}) where T
    return operate!(+, T, operate(-, T, g), f)
end
# Vector Affine +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
    f::MOI.VectorAffineFunction{T}) where T
    return MOI.VectorAffineFunction(operate_terms(op, f.terms),
                      op.(f.constants))
end
function operate(op::Union{typeof(-)}, ::Type{T},
    f::MOI.VectorQuadraticFunction{T}) where T
    return MOI.VectorQuadraticFunction(
        operate_terms(op, f.affine_terms),
        operate_terms(op, f.quadratic_terms),
                      op.(f.constants))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::MOI.VectorAffineFunction{T},
    g::VectorAffineLike{T}) where T
    return operate!(op, T, copy(f), g)
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::MOI.VectorAffineFunction{T},
    g::MOI.VectorQuadraticFunction{T}) where T
    MOI.VectorQuadraticFunction([f.terms; operate_terms(op, g.affine_terms)],
                  operate_terms(op, g.quadratic_terms),
                  op.(f.constants, g.constants))
end

# Vector Quadratic +/- ...
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::MOI.VectorQuadraticFunction{T},
    g::VectorQuadraticLike{T}) where T
    operate!(op, T, copy(f), g)
end

function Base.:+(args::VectorLike{T}...) where T
    return operate(+, T, args...)
end
# Base.:+(α::Vector{T}, f::VectorLike{T}...) is too general as it also covers
# Base.:+(α::Vector) which is type piracy
function Base.:+(α::Vector{T}, f::VectorLike{T}, g::VectorLike{T}...) where T
    return operate(+, T, α, f, g...)
end
function Base.:+(f::VectorLike{T}, α::Vector{T}) where T
    return operate(+, T, f, α)
end
function Base.:-(args::VectorLike{T}...) where T
    return operate(-, T, args...)
end
function Base.:-(f::VectorLike{T}, α::Vector{T}) where T
    return operate(-, T, f, α)
end
function Base.:-(α::Vector{T}, f::VectorLike{T}) where T
    return operate(-, T, α, f)
end

####################################### * ######################################
function promote_operation(::typeof(*), ::Type{T}, ::Type{T},
                           ::Type{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}}) where T
    return MOI.ScalarAffineFunction{T}
end
function promote_operation(::typeof(*), ::Type{T},
                           ::Type{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}},
                           ::Type{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}}) where T
    return MOI.ScalarQuadraticFunction{T}
end


function operate!(::typeof(*), ::Type{T}, f::MOI.SingleVariable, α::T) where T
    return operate(*, T, α, f)
end
function operate(::typeof(*), ::Type{T}, α::T, f::MOI.SingleVariable) where T
    MOI.ScalarAffineFunction{T}([MOI.ScalarAffineTerm(α, f.variable)], zero(T))
end

function operate!(::typeof(*), ::Type{T},
                  f::Union{MOI.ScalarAffineFunction{T},
                           MOI.ScalarQuadraticFunction{T}}, α::T) where T
    map_terms!(term -> operate_term(*, α, term), f)
    f.constant *= α
    return f
end
function operate(::typeof(*), ::Type{T}, α::T, f::MOI.ScalarAffineFunction) where T
    return operate!(*, T, copy(f), α)
end

function operate(::typeof(*), ::Type{T}, α::T,
                 f::MOI.ScalarQuadraticFunction) where T
    return operate!(*, T, copy(f), α)
end

function operate(::typeof(*), ::Type{T}, f::MOI.SingleVariable,
                 g::MOI.SingleVariable) where T
    return MOI.ScalarQuadraticFunction(MOI.ScalarAffineTerm{T}[],
                                       [MOI.ScalarQuadraticTerm(one(T),
                                                                f.variable,
                                                                g.variable)],
                                       zero(T))
end

function operate(::typeof(*), ::Type{T}, f::MOI.ScalarAffineFunction{T},
                 g::MOI.SingleVariable) where T
    aff_terms = [MOI.ScalarAffineTerm(f.constant, g.variable)]
    quad_terms = map(t -> MOI.ScalarQuadraticTerm(t.coefficient,
                                                  t.variable_index,
                                                  g.variable),
                     f.terms)
    return MOI.ScalarQuadraticFunction(aff_terms, quad_terms, zero(T))
end

function operate(::typeof(*), ::Type{T}, f::MOI.ScalarAffineFunction{T},
                 g::MOI.ScalarAffineFunction{T}) where T
    nfterms = length(f.terms)
    ngterms = length(g.terms)
    quad_terms = Vector{MOI.ScalarQuadraticTerm{T}}(undef, nfterms * ngterms)
    k = 0
    for t1 in f.terms
        for t2 in g.terms
            k += 1
            quad_terms[k] = operate_term(*, t1, t2)
        end
    end
    @assert k == length(quad_terms)
    if iszero(f.constant)
        if iszero(g.constant)
            aff_terms = MOI.ScalarAffineTerm{T}[]
        else
            aff_terms = operate_term.(*, g.constant, f.terms)
        end
    else
        if iszero(g.constant)
            aff_terms = operate_term.(*, f.constant, g.terms)
        else
            aff_terms = Vector{MOI.ScalarAffineTerm{T}}(undef,
                                                        nfterms + ngterms)
            map!(t -> operate_term(*, g.constant, t), aff_terms, f.terms)
            for i in 1:ngterms
                aff_terms[nfterms + i] = operate_term(*, f.constant, g.terms[i])
            end
        end
    end
    constant = f.constant * g.constant
    return MOI.ScalarQuadraticFunction(aff_terms, quad_terms, constant)
end

function Base.:*(args::ScalarLike{T}...) where T
    return operate(*, T, args...)
end
function Base.:*(f::T, g::ScalarLike{T}) where T
    return operate(*, T, f, g)
end
function Base.:*(f::ScalarLike{T}, g::T) where T
    return operate(*, T, g, f)
end

####################################### / ######################################
function promote_operation(::typeof(/), ::Type{T},
                           ::Type{<:Union{MOI.SingleVariable,
                                          MOI.ScalarAffineFunction{T}}},
                           ::Type{T}) where T
    MOI.ScalarAffineFunction{T}
end
function promote_operation(::typeof(/), ::Type{T},
                           ::Type{MOI.ScalarQuadraticFunction{T}},
                           ::Type{T}) where T
    MOI.ScalarQuadraticFunction{T}
end

function operate!(::typeof(/), ::Type{T}, f::MOI.SingleVariable,
                  α::T) where T
    return operate(/, T, f, α)
end
function operate(::typeof(/), ::Type{T}, f::MOI.SingleVariable,
                 α::T) where T
    return MOI.ScalarAffineFunction{T}([MOI.ScalarAffineTerm(inv(α),
                                                             f.variable)],
                                       zero(T))
end

function operate!(::typeof(/), ::Type{T}, f::MOI.ScalarAffineFunction{T},
                  α::T) where T
    f.terms .= operate_term.(/, f.terms, α)
    f.constant /= α
    return f
end

function operate!(::typeof(/), ::Type{T}, f::MOI.ScalarQuadraticFunction{T},
                  α::T) where T
    f.affine_terms .= operate_term.(/, f.affine_terms, α)
    f.quadratic_terms .= operate_term.(/, f.quadratic_terms, α)
    f.constant /= α
    return f
end

function operate(::typeof(/), ::Type{T},
                 f::Union{MOI.ScalarAffineFunction{T},
                          MOI.ScalarQuadraticFunction{T}}, α::T) where T
    return operate!(/, T, copy(f), α)
end

function Base.:/(args::ScalarLike{T}...) where T
    return operate(/, T, args...)
end
function Base.:/(f::ScalarLike{T}, g::T) where T
    return operate(/, T, f, g)
end

## sum
function operate(::typeof(sum), ::Type{T}, vis::Vector{MOI.VariableIndex}) where T
    return MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(one(T), vis), zero(T))
end

#################### Concatenation of MOI functions: `vcat` ####################
"""
    fill_vector(vector::Vector, ::Type{T}, fill_func::Function,
                dim_func::Function, funcs) where T

Fill the vector `vector` with
`fill_func(vector, vector_offset, output_offset, func)` for each function `func`
in `funcs` where `vector_offset` (resp. `output_offset`) is the sum of
`dim_func(T, func)` (resp. `output_dim(T, func)`) of previous functions of
`func`.

    fill_vector(vector::Vector, ::Type{T}, vector_offset::Int,
                     output_offset::Int, fill_func::Function,
                     dim_func::Function, funcs...) where T

Same than previous method but starting with possible nonzero `vector_offset` and
`output_offset`.
"""
function fill_vector end

function fill_vector(vector::Vector, ::Type{T}, fill_func::Function,
                     dim_func::Function, funcs) where T
    vector_offset = 0
    output_offset = 0
    for func in funcs
        fill_func(vector, vector_offset, output_offset, func)
        vector_offset += dim_func(T, func)
        output_offset += output_dim(T, func)
    end
    @assert length(vector) == vector_offset
end
function fill_vector(vector::Vector, ::Type, vector_offset::Int,
                     output_offset::Int, fill_func::Function,
                     dim_func::Function)
    @assert length(vector) == vector_offset
end
function fill_vector(vector::Vector, ::Type{T}, vector_offset::Int,
                     output_offset::Int, fill_func::Function,
                     dim_func::Function, func, funcs...) where T
    fill_func(vector, vector_offset, output_offset, func)
    fill_vector(vector, T, vector_offset + dim_func(T, func),
                output_offset + output_dim(T, func), fill_func, dim_func,
                funcs...)
end

number_of_affine_terms(::Type{T}, ::T) where T = 0
number_of_affine_terms(::Type, ::SVF) = 1
number_of_affine_terms(::Type, f::VVF) = length(f.variables)
function number_of_affine_terms(
    ::Type{T}, f::Union{SAF{T}, VAF{T}}) where T
    return length(f.terms)
end
function number_of_affine_terms(
    ::Type{T}, f::Union{SQF{T}, VQF{T}}) where T
    return length(f.affine_terms)
end

function number_of_quadratic_terms(
    ::Type{T}, f::Union{SQF{T}, VQF{T}}) where T
    return length(f.quadratic_terms)
end

function offset_term(t::MOI.ScalarAffineTerm, offset::Int)
    return MOI.VectorAffineTerm(offset + 1, t)
end
function offset_term(t::MOI.VectorAffineTerm, offset::Int)
    return MOI.VectorAffineTerm(offset + t.output_index, t.scalar_term)
end
function offset_term(t::MOI.ScalarQuadraticTerm, offset::Int)
    return MOI.VectorQuadraticTerm(offset + 1, t)
end
function offset_term(t::MOI.VectorQuadraticTerm, offset::Int)
    return MOI.VectorQuadraticTerm(offset + t.output_index, t.scalar_term)
end

function fill_terms(terms::Vector{MOI.VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::T) where T
end
function fill_terms(terms::Vector{MOI.VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::SVF) where T
    terms[offset + 1] = offset_term(
        MOI.ScalarAffineTerm(one(T), func.variable), output_offset)
end
function fill_terms(terms::Vector{MOI.VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::VVF) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= MOI.VectorAffineTerm.(
        output_offset .+ (1:n), MOI.ScalarAffineTerm.(one(T), func.variables))
end
function fill_terms(terms::Vector{MOI.VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{SAF{T}, VAF{T}}) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.terms, output_offset)
end
function fill_terms(terms::Vector{MOI.VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{SQF{T}, VQF{T}}) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.affine_terms, output_offset)
end
function fill_terms(terms::Vector{MOI.VectorQuadraticTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{SQF{T}, VQF{T}}) where T
    n = number_of_quadratic_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.quadratic_terms, output_offset)
end

output_dim(::Type{T}, ::T) where T = 1
output_dim(::Type, func::MOI.AbstractFunction) = MOI.output_dimension(func)
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::T) where T
    constant[offset + 1] = func
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{SVF, VVF}) where T
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{SAF{T}, SQF{T}}) where T
    constant[offset + 1] = func.constant
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{VAF{T}, VQF{T}}) where T
    n = MOI.output_dimension(func)
    constant[offset .+ (1:n)] .= func.constants
end

"""
    vectorize(funcs::AbstractVector{MOI.SingleVariable})

Returns the vector of scalar affine functions in the form of a
`MOI.VectorAffineFunction{T}`.
"""
function vectorize(funcs::AbstractVector{MOI.SingleVariable})
    vars = MOI.VariableIndex[func.variable for func in funcs]
    return MOI.VectorOfVariables(vars)
end

"""
    vectorize(funcs::AbstractVector{MOI.ScalarAffineFunction{T}}) where T

Returns the vector of scalar affine functions in the form of a
`MOI.VectorAffineFunction{T}`.
"""
function vectorize(funcs::AbstractVector{MOI.ScalarAffineFunction{T}}) where T
    nterms = sum(func -> number_of_affine_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    terms = Vector{MOI.VectorAffineTerm{T}}(undef, nterms)
    constant = zeros(T, out_dim)
    fill_vector(terms, T, fill_terms, number_of_affine_terms, funcs)
    fill_vector(constant, T, fill_constant, output_dim, funcs)
    return VAF(terms, constant)
end

"""
    vectorize(funcs::AbstractVector{MOI.ScalarQuadraticFunction{T}}) where T

Returns the vector of scalar quadratic functions in the form of a
`MOI.VectorQuadraticFunction{T}`.
"""
function vectorize(funcs::AbstractVector{MOI.ScalarQuadraticFunction{T}}) where T
    num_affine_terms = sum(func -> number_of_affine_terms(T, func), funcs)
    num_quadratic_terms = sum(func -> number_of_quadratic_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    affine_terms = Vector{MOI.VectorAffineTerm{T}}(undef, num_affine_terms)
    quadratic_terms = Vector{MOI.VectorQuadraticTerm{T}}(undef, num_quadratic_terms)
    constant = zeros(T, out_dim)
    fill_vector(affine_terms, T, fill_terms, number_of_affine_terms, funcs)
    fill_vector(quadratic_terms, T, fill_terms, number_of_quadratic_terms, funcs)
    fill_vector(constant, T, fill_constant, output_dim, funcs)
    return VQF(affine_terms, quadratic_terms, constant)
end

function promote_operation(::typeof(vcat), ::Type{T},
                           ::Type{<:Union{ScalarAffineLike{T}, VVF, VAF{T}}}...) where T
    return VAF{T}
end
function promote_operation(
    ::typeof(vcat), ::Type{T},
    ::Type{<:Union{ScalarQuadraticLike{T}, VVF, VAF{T}, VQF{T}}}...) where T
    return VQF{T}
end

function operate(::typeof(vcat), ::Type{T},
                 funcs::Union{ScalarAffineLike{T}, VVF, VAF{T}}...) where T
    nterms = sum(func -> number_of_affine_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    terms = Vector{MOI.VectorAffineTerm{T}}(undef, nterms)
    constant = zeros(T, out_dim)
    fill_vector(terms, T, 0, 0, fill_terms, number_of_affine_terms, funcs...)
    fill_vector(constant, T, 0, 0, fill_constant, output_dim, funcs...)
    return VAF(terms, constant)
end
function operate(::typeof(vcat), ::Type{T},
                 funcs::Union{ScalarQuadraticLike{T}, VVF, VAF{T}, VQF{T}}...) where T
    num_affine_terms = sum(func -> number_of_affine_terms(T, func), funcs)
    num_quadratic_terms = sum(func -> number_of_quadratic_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    affine_terms = Vector{MOI.VectorAffineTerm{T}}(undef, num_affine_terms)
    quadratic_terms = Vector{MOI.VectorQuadraticTerm{T}}(undef, num_quadratic_terms)
    constant = zeros(T, out_dim)
    fill_vector(affine_terms, T, 0, 0, fill_terms, number_of_affine_terms, funcs...)
    fill_vector(quadratic_terms, T, 0, 0, fill_terms, number_of_quadratic_terms, funcs...)
    fill_vector(constant, T, 0, 0, fill_constant, output_dim, funcs...)
    return VQF(affine_terms, quadratic_terms, constant)
end
