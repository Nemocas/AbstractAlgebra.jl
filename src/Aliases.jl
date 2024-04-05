# alternative names for some functions from Base
@alias is_empty isempty
@alias is_equal isequal
@alias is_even iseven
@alias is_finite isfinite
@alias is_inf isinf
@alias is_integer isinteger
@alias is_less isless
@alias is_odd isodd
@alias is_one isone
@alias is_real isreal
@alias is_subset issubset
@alias is_valid isvalid
@alias is_zero iszero
@alias number_of_digits ndigits

# TODO: next breaking release: remove the if guard around the block
if @__MODULE__() == AbstractAlgebra

    # alternative names for some functions from LinearAlgebra
    # we don't use the `@alias` macro here because we provide custom
    # docstrings for these aliases
    const is_diagonal = isdiag
    const is_hermitian = ishermitian
    const is_symmetric = issymmetric
    const is_lower_triangular = istril
    const is_upper_triangular = istriu
end


# predeclare some functions to allow defining aliases for some of our own functions
function number_of_columns end
function number_of_generators end
function number_of_rows end
function number_of_variables end

@alias ncols number_of_columns
@alias ngens number_of_generators
@alias nrows number_of_rows
@alias nvars number_of_variables
