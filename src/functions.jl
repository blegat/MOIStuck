abstract type AbstractFunction end
function output_dimension end
abstract type AbstractScalarFunction <: AbstractFunction end
output_dimension(::AbstractScalarFunction) = 1
abstract type AbstractVectorFunction <: AbstractFunction end
struct SingleVariable <: AbstractScalarFunction
    variable::VariableIndex
end
struct VectorOfVariables <: AbstractVectorFunction
    variables::Vector{VariableIndex}
end
output_dimension(f::VectorOfVariables) = length(f.variables)
struct ScalarAffineTerm{T}
    coefficient::T
    variable_index::VariableIndex
end
mutable struct ScalarAffineFunction{T} <: AbstractScalarFunction
    terms::Vector{ScalarAffineTerm{T}}
    constant::T
end
struct VectorAffineTerm{T}
    output_index::Int64
    scalar_term::ScalarAffineTerm{T}
end
function VectorAffineTerm(output_index::Base.Integer, scalar_term::ScalarAffineTerm)
    VectorAffineTerm(convert(Int64, output_index), scalar_term)
end
struct VectorAffineFunction{T} <: AbstractVectorFunction
    terms::Vector{VectorAffineTerm{T}}
    constants::Vector{T}
end
output_dimension(f::VectorAffineFunction) = length(f.constants)
struct ScalarQuadraticTerm{T}
    coefficient::T
    variable_index_1::VariableIndex
    variable_index_2::VariableIndex
end
mutable struct ScalarQuadraticFunction{T} <: AbstractScalarFunction
    affine_terms::Vector{ScalarAffineTerm{T}}
    quadratic_terms::Vector{ScalarQuadraticTerm{T}}
    constant::T
end
struct VectorQuadraticTerm{T}
    output_index::Int64
    scalar_term::ScalarQuadraticTerm{T}
end
function VectorQuadraticTerm(output_index::Base.Integer, scalar_term::ScalarQuadraticTerm)
    VectorQuadraticTerm(convert(Int64, output_index), scalar_term)
end


"""
    VectorQuadraticFunction{T}(affine_terms, quadratic_terms, constant)

The vector-valued quadratic function with i`th` component ("output index")
defined as ``\\frac{1}{2}x^TQ_ix + a_i^T x + b_i``, where:
* ``a_i`` is a sparse vector specified by the `VectorAffineTerm`s with
  `output_index == i`.
* ``b_i`` is a scalar specified by `constants[i]`
* ``Q_i`` is a symmetric matrix specified by the `VectorQuadraticTerm` with
  `output_index == i`.

Duplicate indices in ``a_i`` or ``Q_i`` are accepted, and the corresponding
coefficients are summed together. "Mirrored" indices `(q,r)` and `(r,q)` (where
`r` and `q` are `VariableIndex`es) are considered duplicates; only one need be
specified.
"""
struct VectorQuadraticFunction{T} <: AbstractVectorFunction
    affine_terms::Vector{VectorAffineTerm{T}}
    quadratic_terms::Vector{VectorQuadraticTerm{T}}
    constants::Vector{T}
end
output_dimension(f::VectorQuadraticFunction) = length(f.constants)

# Function modifications


"""
    AbstractFunctionModification

An abstract supertype for structs which specify partial modifications to functions, to be used for making small modifications instead of replacing the functions entirely.
"""
abstract type AbstractFunctionModification end

"""
    ScalarConstantChange{T}(new_constant::T)

A struct used to request a change in the constant term of a scalar-valued function.
Applicable to `ScalarAffineFunction` and `ScalarQuadraticFunction`.
"""
struct ScalarConstantChange{T} <: AbstractFunctionModification
    new_constant::T
end

"""
    VectorConstantChange{T}(new_constant::Vector{T})

A struct used to request a change in the constant vector of a vector-valued function.
Applicable to `VectorAffineFunction` and `VectorQuadraticFunction`.
"""
struct VectorConstantChange{T} <: AbstractFunctionModification
    new_constant::Vector{T}
end

"""
    ScalarCoefficientChange{T}(variable::VariableIndex, new_coefficient::T)

A struct used to request a change in the linear coefficient of a single variable
in a scalar-valued function.
Applicable to `ScalarAffineFunction` and `ScalarQuadraticFunction`.
"""
struct ScalarCoefficientChange{T} <: AbstractFunctionModification
    variable::VariableIndex
    new_coefficient::T
end

# Note: MultiRowChange is mutable because its `variable` field of an immutable
# type, while `new_coefficients` is of a mutable type, meaning that creating a `MultiRowChange`
# allocates, and it is desirable to provide a zero-allocation option for working with
# MultiRowChanges. See https://github.com/JuliaOpt/MathOptInterface.jl/pull/343.
"""
    MultirowChange{T}(variable::VariableIndex, new_coefficients::Vector{Tuple{Int64, T}})

A struct used to request a change in the linear coefficients of a single
variable in a vector-valued function. New coefficients are specified by
`(output_index, coefficient)` tuples. Applicable to `VectorAffineFunction` and
`VectorQuadraticFunction`.
"""
mutable struct MultirowChange{T} <: AbstractFunctionModification
    variable::VariableIndex
    new_coefficients::Vector{Tuple{Int64, T}}
end

function MultirowChange(variable::VariableIndex, new_coefficients::Vector{Tuple{Ti, T}}) where {Ti<:Base.Integer, T}
    MultirowChange(variable, [(convert(Int64, i), j) for (i,j) in new_coefficients])
end

# Implementation of comparison for functions
Base.:(==)(f::VectorOfVariables, g::VectorOfVariables) = f.variables == g.variables

Base.isapprox(f::Union{SingleVariable, VectorOfVariables}, g::Union{SingleVariable, VectorOfVariables}; kwargs...) = f == g

# For affine and quadratic functions, terms are compressed in a dictionary using `_dicts` and then the dictionaries are compared with `dict_isapprox`
function dict_isapprox(d1::Dict, d2::Dict{<:Any, T}; kwargs...) where T
    all(kv -> isapprox(kv.second, Base.get(d2, kv.first, zero(T)); kwargs...), d1)
end

# Build a dictionary where the duplicate keys are summed
function sum_dict(kvs::Vector{Pair{K, V}}) where {K, V}
    d = Dict{K, V}()
    for kv in kvs
        key = kv.first
        d[key] = kv.second + Base.get(d, key, zero(V))
    end
    d
end

"""
    coefficient(t::Union{ScalarAffineTerm, ScalarQuadraticTerm
                         VectorAffineTerm, VectorQuadraticTerm})

Finds the coefficient stored in the term `t`.
"""
function coefficient end

function coefficient(t::Union{ScalarAffineTerm, ScalarQuadraticTerm})
    return t.coefficient
end
function coefficient(t::Union{VectorAffineTerm, VectorQuadraticTerm})
    return t.scalar_term.coefficient
end

"""
    term_indices(t::Union{ScalarAffineTerm, ScalarQuadraticTerm,
                         VectorAffineTerm, VectorQuadraticTerm})

Returns the indices of the input term `t` as a tuple of `Int`s.

* For `t::ScalarAffineTerm`, this is a 1-tuple of the variable index.
* For `t::ScalarQuadraticTerm`, this is a 2-tuple of the variable indices
  in non-decreasing order.
* For `t::VectorAffineTerm`, this is a 2-tuple of the row/output and
  variable indices.
* For `t::VectorQuadraticTerm`, this is a 3-tuple of the row/output and
  variable indices in non-decreasing order.
"""
term_indices(t::ScalarAffineTerm) = (t.variable_index.value,)
function term_indices(t::ScalarQuadraticTerm)
    return minmax(t.variable_index_1.value, t.variable_index_2.value)
end
function term_indices(t::Union{VectorAffineTerm, VectorQuadraticTerm})
    return (t.output_index, term_indices(t.scalar_term)...)
end

"""
    term_pair(t::Union{ScalarAffineTerm, ScalarQuadraticTerm,
                       VectorAffineTerm, VectorQuadraticTerm})

Returns the pair [`term_indices`](@ref) `=>` [`coefficient`](@ref) of the term.
"""
function term_pair(t::Union{ScalarAffineTerm, ScalarQuadraticTerm,
                            VectorAffineTerm, VectorQuadraticTerm})
    term_indices(t) => coefficient(t)
end

function _dicts(f::Union{ScalarAffineFunction, VectorAffineFunction})
    return (sum_dict(term_pair.(f.terms)),)
end

function _dicts(f::Union{ScalarQuadraticFunction, VectorQuadraticFunction})
    return (sum_dict(term_pair.(f.affine_terms)),
            sum_dict(term_pair.(f.quadratic_terms)))
end

"""
    constant(f::Union{ScalarAffineFunction, ScalarQuadraticFunction})

Returns the constant term of the scalar function
"""
constant(f::Union{ScalarAffineFunction, ScalarQuadraticFunction}) = f.constant

"""
    constant(f::Union{VectorAffineFunction, VectorQuadraticFunction})

Returns the vector of constant terms of the vector function
"""
constant(f::Union{VectorAffineFunction, VectorQuadraticFunction}) = f.constants

function Base.isapprox(f::F, g::G; kwargs...) where {F<:Union{ScalarAffineFunction, ScalarQuadraticFunction, VectorAffineFunction, VectorQuadraticFunction},
                                                     G<:Union{ScalarAffineFunction, ScalarQuadraticFunction, VectorAffineFunction, VectorQuadraticFunction}}
    isapprox(constant(f), constant(g); kwargs...) && all(dict_isapprox.(_dicts(f), _dicts(g); kwargs...))
end

constant(f::Union{ScalarAffineFunction, ScalarQuadraticFunction}, T::DataType) = constant(f)
constant(f::Union{VectorAffineFunction, VectorQuadraticFunction}, T::DataType) = constant(f)

"""
    constant(f::SingleVariable, T::DataType)

The constant term of a `SingleVariable` function is
the zero value of the specified type `T`.
"""
constant(f::SingleVariable, T::DataType) = zero(T)

"""
    constant(f::VectorOfVariables, T::DataType)

The constant term of a `VectorOfVariables` function is a
vector of zero values of the specified type `T`.
"""
constant(f::VectorOfVariables, T::DataType) = zeros(T, length(f.variables))

# isbits type, nothing to copy
Base.copy(func::SingleVariable) = func

Base.copy(func::VectorOfVariables) = VectorOfVariables(copy(func.variables))

"""
    copy(func::Union{ScalarAffineFunction, VectorAffineFunction})

Return a new affine function with a shallow copy of the terms and constant(s)
from `func`.
"""
Base.copy(func::F) where {F <: Union{ScalarAffineFunction, VectorAffineFunction}} = F(copy(func.terms), copy(constant(func)))

"""
    copy(func::Union{ScalarQuadraticFunction, VectorQuadraticFunction})

Return a new quadratic function with a shallow copy of the terms and constant(s)
from `func`.
"""
Base.copy(func::F) where {F <: Union{ScalarQuadraticFunction, VectorQuadraticFunction}} = F(copy(func.affine_terms), copy(func.quadratic_terms), copy(constant(func)))

function ScalarAffineFunction{T}(f::SingleVariable) where T
    ScalarAffineFunction([ScalarAffineTerm(one(T), f.variable)], zero(T))
end
function VectorAffineFunction{T}(f::VectorOfVariables) where T
    n = length(f.variables)
    return VectorAffineFunction(map(i -> VectorAffineTerm(i, ScalarAffineTerm(one(T), f.variables[i])), 1:n), zeros(T, n))
end
function Base.convert(::Type{SingleVariable}, f::ScalarAffineFunction)
    if !iszero(f.constant) || !isone(length(f.terms)) || !isone(f.terms[1].coefficient)
        throw(InexactError(:convert, SingleVariable, f))
    end
    return SingleVariable(f.terms[1].variable_index)
end
function Base.convert(::Type{SingleVariable},
                      f::ScalarQuadraticFunction{T}) where T
    return convert(SingleVariable, convert(ScalarAffineFunction{T}, f))
end
