import sympy
import spb
from sympy.matrices.expressions import CompanionMatrix
from sympy.integrals import integrate

## GLOBALS
BASIS_ID_STRINGS = { "Bernstein": "B",
                     "Chebyshev": "T",
                     "Lagrange":  "L",
                     "Legendre":  "P",
                     "Monomial":  "M" }

def ChangeOfVariable( from_value, from_domain, to_domain ):
    ## SHIFT TO ZERO
    to_value = from_value - from_domain[0]
    ## SCALE TO EQUAL RANGE
    to_value *= ( to_domain[1] - to_domain[0] ) / ( from_domain[1] - from_domain[0] )
    ## APPLY TO DOMAIN
    to_value += to_domain[0]
    return to_value

def PolynomialRoots( polynomial: sympy.Poly ):
    polynomial_coeffs = polynomial.all_coeffs()
    monic_polynomial = sympy.poly( polynomial / polynomial_coeffs[0] )
    companion_matrix = sympy.Matrix( CompanionMatrix( monic_polynomial ) )
    eigenvalues = sorted( list( companion_matrix.eigenvals().keys() ) )
    return eigenvalues

## POLYNOMIAL BASIS FUNCTIONS
def MonomialBasis( degree, variate, domain ):
    variate = ChangeOfVariable( variate, domain, [0, 1] )
    basis = []
    for n in range( 0, degree + 1 ):
        basis.append( sympy.poly( variate ** ( n ), variate ) )
    return basis

def LagrangeBasis( degree, variate, domain ):
    step = ( domain[1] - domain[0] ) / ( degree )
    node = [ domain[0] + i*step for i in range( degree + 1 ) ]
    basis = []
    for i in range( 0, degree + 1 ):
        basis.append( sympy.poly( variate ** 0, variate ) )
        for j in range( 0, degree + 1 ):
            if i != j:
                basis[i] *= ( ( variate - node[j] ) / ( node[i] - node[j] ) )
    return basis

def BernsteinBasis( degree, variate, domain ):
    variate = ChangeOfVariable( variate, domain, [0, 1] )
    basis = []
    for i in range( 0, degree + 1 ):
        bin_coeff = sympy.binomial( degree, i )
        monomial = sympy.poly( variate ** i, variate )
        polynomial = sympy.poly( ( ( 1 - variate ) ** ( degree - i ) ), variate )
        basis.append(  bin_coeff * monomial * polynomial )
    return basis

def LegendreBasis( degree, variate, domain ):
    variate = ChangeOfVariable( variate, domain, [-1, 1] )
    basis = []
    for i in range( 0, degree + 1 ):
        if i == 0:
            basis.append( sympy.poly( variate ** 0, variate ) )
        elif i == 1:
            basis.append( sympy.poly( variate ** 1, variate ) )
        else:
            n = i - 1
            term_1 = ( 2 * n + 1 ) * variate * basis[n]
            term_2 = n * basis[n-1]
            basis.append( ( term_1 - term_2 ) / ( n + 1 ) )
    return basis

def ChebyshevBasis( degree, variate, domain ):
    variate = ChangeOfVariable( variate, domain, [-1, 1] )
    basis = []
    for i in range( 0, degree + 1 ):
        if i == 0:
            basis.append( sympy.poly( variate ** 0, variate ) )
        elif i == 1:
            basis.append( sympy.poly( variate ** 1, variate ) )
        else:
            basis.append( ( 2 * variate * basis[i-1] ) - basis[i-2] )
    return basis

def PolynomialBasisFunction( basis_name, degree, variate, domain ):
    match basis_name:
        case "Bernstein":
            basis = BernsteinBasis( degree, variate, domain )
        case "Chebyshev":
            basis = ChebyshevBasis( degree, variate, domain )
        case "Lagrange":
            basis = LagrangeBasis( degree, variate, domain )
        case "Legendre":
            basis = LegendreBasis( degree, variate, domain )
        case "Monomial":
            basis = MonomialBasis( degree, variate, domain )
    return basis

## Change of Basis
def VectorChangeOfBasis( from_basis, to_basis, from_coeffs ):
    basis_dimension = len( from_coeffs )
    D = sympy.zeros( basis_dimension, basis_dimension )
    C = sympy.zeros( basis_dimension, basis_dimension )
    for i in range( 0, basis_dimension ):
        for j in range( 0, basis_dimension ):
            D[i,j] = to_basis[:,i].dot( to_basis[:, j] )
            C[i,j] = to_basis[:,i].dot( from_basis[:,j] )
    R = D.solve( C )
    to_coeffs = R.multiply( from_coeffs )
    return to_coeffs, R

def PolynomialChangeOfBasis( from_basis, to_basis, from_coeffs, domain ):
    basis_dim = len( from_coeffs )
    D = sympy.zeros( basis_dim )
    C = sympy.zeros( basis_dim )
    variate = list( from_basis[0].atoms( sympy.Symbol ) )[0]
    for i in range( 0, basis_dim ):
        for j in range( 0, basis_dim ):
            integrand = to_basis[i].as_expr() * to_basis[j].as_expr()
            D[i,j] = integrate( integrand, ( variate, domain[0], domain[1] ) )
            integrand = to_basis[i].as_expr() * from_basis[j].as_expr()
            C[i,j] = integrate( integrand, ( variate, domain[0], domain[1] ) )
    R = D.solve( C )
    to_coeffs = R.multiply( from_coeffs )

## Plotting
def PlotPolynomialBasis( basis_name, degree, variate, domain ):
    basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
    basis_id_string = BASIS_ID_STRINGS[basis_name]
    for n in range( 0, len( basis ) ):
        if n == 0:
            plt = spb.plot( basis[n].as_expr(), (variate, domain[0], domain[1] ), label=f"${basis_id_string}_{n}(x)$", show=False )
        else:
            plt += spb.plot( basis[n].as_expr(), (variate, domain[0], domain[1] ), label=f"${basis_id_string}_{n}(x)$", show=False )
    plt.show()

## UTILITIES
def BasisToLatexString( basis_array, basis_id_string ):
    latex_string = basis_id_string + "=" + "\\begin{bmatrix}"
    for basis in basis_array:
        latex_string += f"{sympy.latex( sympy.simplify( basis.as_expr() ) )} \\\\"
    latex_string += "\\end{bmatrix}"
    return latex_string
