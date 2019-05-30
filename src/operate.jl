# Functions convertible to a ScalarAffineFunction
const ScalarAffineLike{T} = Union{T, SingleVariable, ScalarAffineFunction{T}}
# Functions convertible to a ScalarQuadraticFunction
const ScalarQuadraticLike{T} = Union{ScalarAffineLike{T}, ScalarQuadraticFunction{T}}

# Used for overloading Base operator functions so `T` is not in the union to
# avoid overloading e.g. `+(::Float64, ::Float64)`
const ScalarLike{T} = Union{SingleVariable, ScalarAffineFunction{T},
                            ScalarQuadraticFunction{T}}

# Functions convertible to a VectorAffineFunction
const VectorAffineLike{T} = Union{Vector{T}, VectorOfVariables, VectorAffineFunction{T}}
# Functions convertible to a VectorQuadraticFunction
const VectorQuadraticLike{T} = Union{VectorAffineLike{T}, VectorQuadraticFunction{T}}

# Used for overloading Base operator functions so `T` is not in the union to
# avoid overloading e.g. `+(::Float64, ::Float64)`
const VectorLike{T} = Union{VectorOfVariables, VectorAffineFunction{T},
                            VectorQuadraticFunction{T}}

###################################### +/- #####################################
## promote_operation

function promote_operation(::typeof(-), ::Type{T},
                           ::Type{<:ScalarAffineLike{T}}) where T
    return ScalarAffineFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:ScalarAffineLike{T}},
                           ::Type{<:ScalarAffineLike{T}}) where T
    return ScalarAffineFunction{T}
end
function promote_operation(::typeof(-), ::Type{T},
                           ::Type{<:ScalarQuadraticLike{T}}) where T
    return ScalarQuadraticFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:ScalarQuadraticLike{T}},
                           ::Type{<:ScalarQuadraticLike{T}}) where T
    return ScalarQuadraticFunction{T}
end

## operate!
# + with at least 3 arguments
function operate!(op::typeof(+), ::Type{T}, f, g, h, args...) where T
    operate!(op, T, f, g)
    return operate!(+, T, f, h, args...)
end

# Unary -
function operate!(op::typeof(-), ::Type{T}, f::ScalarQuadraticFunction{T}) where T
    operate_terms!(-, f.quadratic_terms)
    operate_terms!(-, f.affine_terms)
    f.constant = -f.constant
    return f
end


# Scalar Variable +/- ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::SingleVariable,
                  g::ScalarQuadraticLike) where T
    return operate(op, T, f, g)
end
# Scalar Affine +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarAffineFunction{T},
                  g::T) where T
    f.constant = op(f.constant, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarAffineFunction{T},
                  g::SingleVariable) where T
    push!(f.terms, ScalarAffineTerm(op(one(T)), g.variable))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarAffineFunction{T},
                  g::ScalarAffineFunction{T}) where T
    append!(f.terms, operate_terms(op, g.terms))
    f.constant = op(f.constant, g.constant)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarAffineFunction{T},
                  g::ScalarQuadraticFunction{T}) where T
    return operate(op, T, f, g)
end
# Scalar Quadratic +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarQuadraticFunction{T},
                  g::T) where T
    f.constant = op(f.constant, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarQuadraticFunction{T},
                  g::SingleVariable) where T
    push!(f.affine_terms, ScalarAffineTerm(op(one(T)), g.variable))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarQuadraticFunction{T},
                  g::ScalarAffineFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.terms))
    f.constant = op(f.constant, g.constant)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::ScalarQuadraticFunction{T},
                  g::ScalarQuadraticFunction{T}) where T
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
function operate(op::typeof(+), ::Type{T}, f::AbstractFunction) where T
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
function operate(::typeof(-), ::Type{T}, f::SingleVariable) where T
    return ScalarAffineFunction{T}(
        [ScalarAffineTerm(-one(T), f.variable)], zero(T))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::SingleVariable, α::T) where T
    return ScalarAffineFunction{T}(
        [ScalarAffineTerm(one(T), f.variable)], op(α))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::SingleVariable,
                 g::SingleVariable) where T
    return ScalarAffineFunction{T}(
        [ScalarAffineTerm(one(T), f.variable),
         ScalarAffineTerm(op(one(T)), g.variable)],
        zero(T))
end
function operate(op::typeof(+), ::Type{T},
                 f::SingleVariable,
                 g::Union{ScalarAffineFunction{T},
                          ScalarQuadraticFunction{T}}) where T
    return operate(op, T, g, f)
end
function operate(op::typeof(-), ::Type{T},
                 f::SingleVariable,
                 g::Union{ScalarAffineFunction{T},
                          ScalarQuadraticFunction{T}}) where T
    return operate!(+, T, operate(-, T, g), f)
end
# Scalar Affine +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
                 f::ScalarAffineFunction{T}) where T
    return ScalarAffineFunction(operate_terms(op, f.terms),
                                    op(f.constant))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::ScalarAffineFunction{T},
                 g::ScalarAffineLike{T}) where T
    return operate!(op, T, copy(f), g)
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::ScalarAffineFunction{T},
                 g::ScalarQuadraticFunction{T}) where T
    ScalarQuadraticFunction([f.terms; operate_terms(op, g.affine_terms)],
                                operate_terms(op, g.quadratic_terms),
                                op(f.constant, g.constant))
end
# Scalar Quadratic +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
                 f::ScalarQuadraticFunction{T}) where T
    return ScalarQuadraticFunction(
        operate_terms(op, f.affine_terms),
        operate_terms(op, f.quadratic_terms),
        op(f.constant))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::ScalarQuadraticFunction{T},
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
    return VectorAffineFunction{T}
end
function promote_operation(::typeof(-), ::Type{T},
    ::Type{<:VectorQuadraticLike{T}}) where T
    return VectorQuadraticFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:VectorAffineLike{T}},
                           ::Type{<:VectorAffineLike{T}}) where T
    return VectorAffineFunction{T}
end
function promote_operation(::Union{typeof(+), typeof(-)}, ::Type{T},
                           ::Type{<:VectorQuadraticLike{T}},
                           ::Type{<:VectorQuadraticLike{T}}) where T
    return VectorQuadraticFunction{T}
end

# Vector Variable +/- ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorOfVariables,
                  g::VectorQuadraticLike) where T
    return operate(op, T, f, g)
end
# Vector Affine +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorAffineFunction{T},
                  g::Vector{T}) where T
    @assert output_dimension(f) == length(g)
    f.constants .= op.(f.constants, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorAffineFunction{T},
                  g::VectorOfVariables) where T
    d = output_dimension(g)
    @assert output_dimension(f) == d
    append!(f.terms, VectorAffineTerm.(
                         collect(1:d),
                         ScalarAffineTerm.(op(one(T)), g.variables)))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorAffineFunction{T},
                  g::VectorAffineFunction{T}) where T
    append!(f.terms, operate_terms(op, g.terms))
    f.constants .= op.(f.constants, g.constants)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorAffineFunction{T},
                  g::VectorQuadraticFunction{T}) where T
    return operate(op, T, f, g)
end
# Vector Quadratic +/-! ...
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorQuadraticFunction{T},
                  g::Vector{T}) where T
    @assert output_dimension(f) == length(g)
    f.constants .= op.(f.constants, g)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorQuadraticFunction{T},
                  g::VectorOfVariables) where T
    d = output_dimension(g)
    @assert output_dimension(f) == d
    append!(f.affine_terms, VectorAffineTerm.(collect(1:d), ScalarAffineTerm.(op(one(T)), g.variables)))
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorQuadraticFunction{T},
                  g::VectorAffineFunction{T}) where T
    append!(f.affine_terms, operate_terms(op, g.terms))
    f.constants .= op.(f.constants, g.constants)
    return f
end
function operate!(op::Union{typeof(+), typeof(-)}, ::Type{T},
                  f::VectorQuadraticFunction{T},
                  g::VectorQuadraticFunction{T}) where T
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

function operate(op::typeof(-), ::Type{T}, f::VectorOfVariables) where T
    d = output_dimension(f)
    return VectorAffineFunction{T}(
               VectorAffineTerm.(
                   collect(1:d),
                   ScalarAffineTerm.(-one(T), f.variables)),
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
    f::VectorOfVariables, α::Vector{T}) where T
    d = output_dimension(f)
    @assert length(α) == d
    return VectorAffineFunction{T}(
               VectorAffineTerm.(
                   collect(1:d),
                   ScalarAffineTerm.(one(T), f.variables)),
               op.(α))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
                 f::VectorOfVariables,
                 g::VectorOfVariables) where T
    d = output_dimension(f)
    @assert output_dimension(g) == d
    return VectorAffineFunction{T}(
               vcat(
                   VectorAffineTerm.(
                       collect(1:d),
                       ScalarAffineTerm.(one(T), f.variables)),
                   VectorAffineTerm.(
                       collect(1:d),
                       ScalarAffineTerm.(op(one(T)), g.variables))),
               fill(zero(T),d))
end
function operate(op::typeof(+), ::Type{T},
    f::VectorOfVariables,
    g::Union{VectorAffineFunction{T},
                VectorQuadraticFunction{T}}) where T
    return operate(op, T, g, f)
end
function operate(op::typeof(-), ::Type{T},
    f::VectorOfVariables,
    g::Union{VectorAffineFunction{T},
                VectorQuadraticFunction{T}}) where T
    return operate!(+, T, operate(-, T, g), f)
end
# Vector Affine +/- ...
function operate(op::Union{typeof(-)}, ::Type{T},
    f::VectorAffineFunction{T}) where T
    return VectorAffineFunction(operate_terms(op, f.terms),
                      op.(f.constants))
end
function operate(op::Union{typeof(-)}, ::Type{T},
    f::VectorQuadraticFunction{T}) where T
    return VectorQuadraticFunction(
        operate_terms(op, f.affine_terms),
        operate_terms(op, f.quadratic_terms),
                      op.(f.constants))
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::VectorAffineFunction{T},
    g::VectorAffineLike{T}) where T
    return operate!(op, T, copy(f), g)
end
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::VectorAffineFunction{T},
    g::VectorQuadraticFunction{T}) where T
    VectorQuadraticFunction([f.terms; operate_terms(op, g.affine_terms)],
                  operate_terms(op, g.quadratic_terms),
                  op.(f.constants, g.constants))
end

# Vector Quadratic +/- ...
function operate(op::Union{typeof(+), typeof(-)}, ::Type{T},
    f::VectorQuadraticFunction{T},
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
                           ::Type{<:Union{SingleVariable,
                                          ScalarAffineFunction{T}}}) where T
    return ScalarAffineFunction{T}
end
function promote_operation(::typeof(*), ::Type{T},
                           ::Type{<:Union{SingleVariable,
                                          ScalarAffineFunction{T}}},
                           ::Type{<:Union{SingleVariable,
                                          ScalarAffineFunction{T}}}) where T
    return ScalarQuadraticFunction{T}
end


function operate!(::typeof(*), ::Type{T}, f::SingleVariable, α::T) where T
    return operate(*, T, α, f)
end
function operate(::typeof(*), ::Type{T}, α::T, f::SingleVariable) where T
    ScalarAffineFunction{T}([ScalarAffineTerm(α, f.variable)], zero(T))
end

function operate!(::typeof(*), ::Type{T},
                  f::Union{ScalarAffineFunction{T},
                           ScalarQuadraticFunction{T}}, α::T) where T
    map_terms!(term -> operate_term(*, α, term), f)
    f.constant *= α
    return f
end
function operate(::typeof(*), ::Type{T}, α::T, f::ScalarAffineFunction) where T
    return operate!(*, T, copy(f), α)
end

function operate(::typeof(*), ::Type{T}, α::T,
                 f::ScalarQuadraticFunction) where T
    return operate!(*, T, copy(f), α)
end

function operate(::typeof(*), ::Type{T}, f::SingleVariable,
                 g::SingleVariable) where T
    return ScalarQuadraticFunction(ScalarAffineTerm{T}[],
                                       [ScalarQuadraticTerm(one(T),
                                                                f.variable,
                                                                g.variable)],
                                       zero(T))
end

function operate(::typeof(*), ::Type{T}, f::ScalarAffineFunction{T},
                 g::SingleVariable) where T
    aff_terms = [ScalarAffineTerm(f.constant, g.variable)]
    quad_terms = map(t -> ScalarQuadraticTerm(t.coefficient,
                                                  t.variable_index,
                                                  g.variable),
                     f.terms)
    return ScalarQuadraticFunction(aff_terms, quad_terms, zero(T))
end

function operate(::typeof(*), ::Type{T}, f::ScalarAffineFunction{T},
                 g::ScalarAffineFunction{T}) where T
    nfterms = length(f.terms)
    ngterms = length(g.terms)
    quad_terms = Vector{ScalarQuadraticTerm{T}}(undef, nfterms * ngterms)
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
            aff_terms = ScalarAffineTerm{T}[]
        else
            aff_terms = operate_term.(*, g.constant, f.terms)
        end
    else
        if iszero(g.constant)
            aff_terms = operate_term.(*, f.constant, g.terms)
        else
            aff_terms = Vector{ScalarAffineTerm{T}}(undef,
                                                        nfterms + ngterms)
            map!(t -> operate_term(*, g.constant, t), aff_terms, f.terms)
            for i in 1:ngterms
                aff_terms[nfterms + i] = operate_term(*, f.constant, g.terms[i])
            end
        end
    end
    constant = f.constant * g.constant
    return ScalarQuadraticFunction(aff_terms, quad_terms, constant)
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
                           ::Type{<:Union{SingleVariable,
                                          ScalarAffineFunction{T}}},
                           ::Type{T}) where T
    ScalarAffineFunction{T}
end
function promote_operation(::typeof(/), ::Type{T},
                           ::Type{ScalarQuadraticFunction{T}},
                           ::Type{T}) where T
    ScalarQuadraticFunction{T}
end

function operate!(::typeof(/), ::Type{T}, f::SingleVariable,
                  α::T) where T
    return operate(/, T, f, α)
end
function operate(::typeof(/), ::Type{T}, f::SingleVariable,
                 α::T) where T
    return ScalarAffineFunction{T}([ScalarAffineTerm(inv(α),
                                                             f.variable)],
                                       zero(T))
end

function operate!(::typeof(/), ::Type{T}, f::ScalarAffineFunction{T},
                  α::T) where T
    f.terms .= operate_term.(/, f.terms, α)
    f.constant /= α
    return f
end

function operate!(::typeof(/), ::Type{T}, f::ScalarQuadraticFunction{T},
                  α::T) where T
    f.affine_terms .= operate_term.(/, f.affine_terms, α)
    f.quadratic_terms .= operate_term.(/, f.quadratic_terms, α)
    f.constant /= α
    return f
end

function operate(::typeof(/), ::Type{T},
                 f::Union{ScalarAffineFunction{T},
                          ScalarQuadraticFunction{T}}, α::T) where T
    return operate!(/, T, copy(f), α)
end

function Base.:/(args::ScalarLike{T}...) where T
    return operate(/, T, args...)
end
function Base.:/(f::ScalarLike{T}, g::T) where T
    return operate(/, T, f, g)
end

## sum
function operate(::typeof(sum), ::Type{T}, vis::Vector{VariableIndex}) where T
    return ScalarAffineFunction(ScalarAffineTerm.(one(T), vis), zero(T))
end

#################### Concatenation of functions: `vcat` ####################
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
number_of_affine_terms(::Type, ::SingleVariable) = 1
number_of_affine_terms(::Type, f::VectorOfVariables) = length(f.variables)
function number_of_affine_terms(
    ::Type{T}, f::Union{ScalarAffineFunction{T}, VectorAffineFunction{T}}) where T
    return length(f.terms)
end
function number_of_affine_terms(
    ::Type{T}, f::Union{ScalarQuadraticFunction{T}, VectorQuadraticFunction{T}}) where T
    return length(f.affine_terms)
end

function number_of_quadratic_terms(
    ::Type{T}, f::Union{ScalarQuadraticFunction{T}, VectorQuadraticFunction{T}}) where T
    return length(f.quadratic_terms)
end

function offset_term(t::ScalarAffineTerm, offset::Int)
    return VectorAffineTerm(offset + 1, t)
end
function offset_term(t::VectorAffineTerm, offset::Int)
    return VectorAffineTerm(offset + t.output_index, t.scalar_term)
end
function offset_term(t::ScalarQuadraticTerm, offset::Int)
    return VectorQuadraticTerm(offset + 1, t)
end
function offset_term(t::VectorQuadraticTerm, offset::Int)
    return VectorQuadraticTerm(offset + t.output_index, t.scalar_term)
end

function fill_terms(terms::Vector{VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::T) where T
end
function fill_terms(terms::Vector{VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::SingleVariable) where T
    terms[offset + 1] = offset_term(
        ScalarAffineTerm(one(T), func.variable), output_offset)
end
function fill_terms(terms::Vector{VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::VectorOfVariables) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= VectorAffineTerm.(
        output_offset .+ (1:n), ScalarAffineTerm.(one(T), func.variables))
end
function fill_terms(terms::Vector{VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{ScalarAffineFunction{T}, VectorAffineFunction{T}}) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.terms, output_offset)
end
function fill_terms(terms::Vector{VectorAffineTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{ScalarQuadraticFunction{T}, VectorQuadraticFunction{T}}) where T
    n = number_of_affine_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.affine_terms, output_offset)
end
function fill_terms(terms::Vector{VectorQuadraticTerm{T}}, offset::Int,
                    output_offset::Int, func::Union{ScalarQuadraticFunction{T}, VectorQuadraticFunction{T}}) where T
    n = number_of_quadratic_terms(T, func)
    terms[offset .+ (1:n)] .= offset_term.(func.quadratic_terms, output_offset)
end

output_dim(::Type{T}, ::T) where T = 1
output_dim(::Type, func::AbstractFunction) = output_dimension(func)
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::T) where T
    constant[offset + 1] = func
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{SingleVariable, VectorOfVariables}) where T
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{ScalarAffineFunction{T}, ScalarQuadraticFunction{T}}) where T
    constant[offset + 1] = func.constant
end
function fill_constant(constant::Vector{T}, offset::Int,
                       output_offset::Int, func::Union{VectorAffineFunction{T}, VectorQuadraticFunction{T}}) where T
    n = output_dimension(func)
    constant[offset .+ (1:n)] .= func.constants
end

"""
    vectorize(funcs::AbstractVector{SingleVariable})

Returns the vector of scalar affine functions in the form of a
`VectorAffineFunction{T}`.
"""
function vectorize(funcs::AbstractVector{SingleVariable})
    vars = VariableIndex[func.variable for func in funcs]
    return VectorOfVariables(vars)
end

"""
    vectorize(funcs::AbstractVector{ScalarAffineFunction{T}}) where T

Returns the vector of scalar affine functions in the form of a
`VectorAffineFunction{T}`.
"""
function vectorize(funcs::AbstractVector{ScalarAffineFunction{T}}) where T
    nterms = sum(func -> number_of_affine_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    terms = Vector{VectorAffineTerm{T}}(undef, nterms)
    constant = zeros(T, out_dim)
    fill_vector(terms, T, fill_terms, number_of_affine_terms, funcs)
    fill_vector(constant, T, fill_constant, output_dim, funcs)
    return VectorAffineFunction(terms, constant)
end

"""
    vectorize(funcs::AbstractVector{ScalarQuadraticFunction{T}}) where T

Returns the vector of scalar quadratic functions in the form of a
`VectorQuadraticFunction{T}`.
"""
function vectorize(funcs::AbstractVector{ScalarQuadraticFunction{T}}) where T
    num_affine_terms = sum(func -> number_of_affine_terms(T, func), funcs)
    num_quadratic_terms = sum(func -> number_of_quadratic_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    affine_terms = Vector{VectorAffineTerm{T}}(undef, num_affine_terms)
    quadratic_terms = Vector{VectorQuadraticTerm{T}}(undef, num_quadratic_terms)
    constant = zeros(T, out_dim)
    fill_vector(affine_terms, T, fill_terms, number_of_affine_terms, funcs)
    fill_vector(quadratic_terms, T, fill_terms, number_of_quadratic_terms, funcs)
    fill_vector(constant, T, fill_constant, output_dim, funcs)
    return VectorQuadraticFunction(affine_terms, quadratic_terms, constant)
end

function promote_operation(::typeof(vcat), ::Type{T},
                           ::Type{<:Union{ScalarAffineLike{T}, VectorOfVariables, VectorAffineFunction{T}}}...) where T
    return VectorAffineFunction{T}
end
function promote_operation(
    ::typeof(vcat), ::Type{T},
    ::Type{<:Union{ScalarQuadraticLike{T}, VectorOfVariables, VectorAffineFunction{T}, VectorQuadraticFunction{T}}}...) where T
    return VectorQuadraticFunction{T}
end

function operate(::typeof(vcat), ::Type{T},
                 funcs::Union{ScalarAffineLike{T}, VectorOfVariables, VectorAffineFunction{T}}...) where T
    nterms = sum(func -> number_of_affine_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    terms = Vector{VectorAffineTerm{T}}(undef, nterms)
    constant = zeros(T, out_dim)
    fill_vector(terms, T, 0, 0, fill_terms, number_of_affine_terms, funcs...)
    fill_vector(constant, T, 0, 0, fill_constant, output_dim, funcs...)
    return VectorAffineFunction(terms, constant)
end
function operate(::typeof(vcat), ::Type{T},
                 funcs::Union{ScalarQuadraticLike{T}, VectorOfVariables, VectorAffineFunction{T}, VectorQuadraticFunction{T}}...) where T
    num_affine_terms = sum(func -> number_of_affine_terms(T, func), funcs)
    num_quadratic_terms = sum(func -> number_of_quadratic_terms(T, func), funcs)
    out_dim = sum(func -> output_dim(T, func), funcs)
    affine_terms = Vector{VectorAffineTerm{T}}(undef, num_affine_terms)
    quadratic_terms = Vector{VectorQuadraticTerm{T}}(undef, num_quadratic_terms)
    constant = zeros(T, out_dim)
    fill_vector(affine_terms, T, 0, 0, fill_terms, number_of_affine_terms, funcs...)
    fill_vector(quadratic_terms, T, 0, 0, fill_terms, number_of_quadratic_terms, funcs...)
    fill_vector(constant, T, 0, 0, fill_constant, output_dim, funcs...)
    return VectorQuadraticFunction(affine_terms, quadratic_terms, constant)
end
