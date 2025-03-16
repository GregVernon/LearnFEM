module lecture_01_utils

# IMPORTS
using LinearAlgebra
using Symbolics

# EXPORTS
export to_numeric
export diagm
export eigvals
export companion_matrix
export get_coeffs
export polynomial_to_monic

# DEFINITIONS

## TO_NUMERIC
function to_numeric( T::Type, A::Vector{Symbolics.BasicSymbolic} )
    B = zeros( T, length( A ) )
    for idx in eachindex( A )
        B[idx] = T( A[idx].arguments[1].val )
    end
    return B
end

function to_numeric( T::Type, A::Vector{Symbolics.Num} )
    B = zeros( T, length( A ) )
    for idx in eachindex( A )
        B[idx] = T( A[idx].val )
    end
    return B
end

function to_numeric( T::Type, A::Vector{<:Real} )
    B = zeros( T, length( A ) )
    for idx in eachindex( A )
        if typeof( A[idx] ) == Symbolics.Num
            B[idx] = T( A[idx].val )
        else
            B[idx] = T( A[idx] )
        end
    end
    return B
end

function to_numeric( T::Type, A::Matrix{Symbolics.BasicSymbolic} )
    B = zeros( T, size( A, 1 ), size( A, 2 ) )
    for idx in eachindex( A )
        B[idx] = T( A[idx].arguments[1].val )
    end
    return B
end

function to_numeric( T::Type, A::Matrix{Symbolics.Num} )
    B = zeros( T, size( A, 1 ), size( A, 2 ) )
    for idx in eachindex( A )
        B[idx] = T( A[idx].val )
    end
    return B
end

function to_numeric( T::Type, A::Matrix{<:Real} )
    B = zeros( T, size( A, 1 ), size( A, 2 ) )
    for idx in eachindex( A )
        if typeof( A[idx] ) == Symbolics.Num
            B[idx] = T( A[idx].val )
        else
            B[idx] = T( A[idx] )
        end
    end
    return B
end

## DIAGM
function diagm( m::Matrix{BigFloat} )
    num_rows, num_cols = size( m )
    v = zeros( BigFloat, num_rows )
    for idx = 1 : num_rows
        v[idx] = m[idx,idx]
    end
    return v
end

## GET COEFFS
function get_coeffs( polynomial )
    variate = Symbolics.get_variables( polynomial )[1]
    coeffs = []
    p = 0
    while Symbolics.iszero( polynomial ) == false
        monomial = variate ^ p
        append!( coeffs, Symbolics.coeff( polynomial, monomial ) )
        polynomial -= coeffs[end] * monomial
        p += 1
    end
    return coeffs
end

## POLYNOMIAL TO MONIC
function polynomial_to_monic( polynomial )
    variate = Symbolics.get_variables( polynomial )[1]
    coeffs = get_coeffs( polynomial )
    coeffs /= coeffs[end]
    monic_polynomial = Symbolics.Num( 0 )
    for i in eachindex( coeffs )
        monic_polynomial += coeffs[i] * variate^(i-1)
    end
    return monic_polynomial
end

## COMPANION MATRIX
function companion_matrix( polynomial )
    monic_polynomial = polynomial_to_monic( polynomial )
    monic_coeffs = get_coeffs( monic_polynomial )
    poly_degree = length( monic_coeffs ) - 1
    comp_matrix = fill( Symbolics.zero( Symbolics.Num ), poly_degree, poly_degree )
    for i = 0 : poly_degree - 1
        row = i + 1
        if i > 0
            comp_matrix[row,row-1] = 1
        end
        comp_matrix[row, end] = -1 * monic_coeffs[i+1]
    end
    return comp_matrix
end

## EIGVALS
function eigvals(A::Matrix{Symbolics.BasicSymbolic}; tol=1e-10, max_iter=1e4)
    A = to_numeric( BigFloat, A )
    A_prev = A
    A_next = A_prev
    for i in 1:max_iter
        Q, R = qr(A_prev)
        A_next = R * Q
        
        ritz_residual = diagm( A_next ) - diagm( A_prev )
        ritz_residual_norm = norm( ritz_residual )
        if ritz_residual_norm < tol
            break
        end
        
        A_prev = A_next
    end
    
    # Extract eigenvalues from the final matrix symbolically (diagonal elements)
    eigenvalues = diagm( A_next )
    return eigenvalues
end

function eigvals(A::Matrix{Symbolics.Num}; tol=1e-10, max_iter=1e4)
    A = to_numeric( BigFloat, A )
    ritz_values = Vector{Any}( [ diagm( A ) ] )
    for i in 1:max_iter
        Q, R = qr( A )
        A = R * Q
        append!( ritz_values, [ diagm( A ) ] )
        ritz_residual = ritz_values[end] - ritz_values[end-1]
        ritz_residual_norm = norm( ritz_residual, Inf )
        if ritz_residual_norm < tol
            break
        end
    end    
    eigenvalues = ritz_values[end]
    return eigenvalues, ritz_values
end

function eigvals(A::Matrix{BigFloat}; tol=1e-10, max_iter=1e4)
    ritz_values = Vector{any}( [ diagm( A ) ] )
    for i in 1:max_iter
        Q, R = qr( A )
        A = R * Q
        append!( ritz_values, [ diagm( A ) ] )
        ritz_residual = ritz_values[end] - ritz_values[end-1]
        ritz_residual_norm = norm( ritz_residual, Inf )
        # # Check convergence symbolically (approximating with numerical evaluation)
        if ritz_residual_norm < tol
            break
        end
    end
    eigenvalues = ritz_values[end]
    return eigenvalues, ritz_values
end

end