import numpy
import scipy
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

def BasisPolynomialListToExprMatrix( basis ):
    return sympy.Matrix( [ basis_fun.as_expr() for basis_fun in basis ] )

def BasisPolynomialListToLambdaArray( basis, variate ):
    return [ sympy.lambdify( variate, bfun.as_expr() ) for bfun in basis ]

## Scalar Projection
class ScalarProjection:
    def __call__(self, target_fun, basis_name, degree, domain ):
        variate = list( target_fun.atoms( sympy.Symbol ) )[0]
        basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
        basis = BasisPolynomialListToExprMatrix( basis )
        M = self.AssembleGramMatrix( basis, domain, variate )
        F = self.AssembleForceVector( target_fun, basis, domain, variate )
        d = M.solve( F )
        u = ( d.T * basis )[0]
        return u, M, F, basis, d
        
    @staticmethod
    def AssembleGramMatrix( basis, domain, variate ):
        basis_dim = len( basis )
        M = sympy.zeros( basis_dim )
        for i in range( 0, basis_dim ):
            for j in range( 0, basis_dim ):
                M[i,j] = integrate( basis[i] * basis[j], ( variate, domain[0], domain[1] ) )
        return M
    
    @staticmethod
    def AssembleForceVector( target_fun, basis, domain, variate ):
        basis_dim = len( basis )
        F = sympy.zeros( rows=basis_dim, cols=1 )
        for i in range( 0, basis_dim ):
            F[i] = integrate( basis[i] * target_fun, ( variate, domain[0], domain[1] ) )
        return F

class ScalarProjectionFast:
    def __call__(self, target_fun, basis_name, degree, domain ):
        variate = list( target_fun.atoms( sympy.Symbol ) )[0]
        target_fun = sympy.lambdify( variate, target_fun )
        basis = PolynomialBasisFunction( basis_name, degree, variate, domain )
        basis = BasisPolynomialListToLambdaArray( basis, variate )
        M = self.AssembleGramMatrix( basis, domain )
        F = self.AssembleForceVector( target_fun, basis, domain )
        d = numpy.linalg.inv( M ) @ F
        u = lambda x: numpy.sum( [ d[i] * basis[i](x) for i in range( 0, len( basis ) ) ] )
        return u, M, F, basis, d
        
    @staticmethod
    def AssembleGramMatrix( basis, domain ):
        basis_dim = len( basis )
        M = numpy.zeros( shape=( basis_dim, basis_dim ) )
        for i in range( 0, basis_dim ):
            for j in range( 0, basis_dim ):
                M[i,j] = scipy.integrate.quad( lambda x: basis[i](x) * basis[j](x), domain[0], domain[1], limit=1000, epsrel=1e-12 )[0]
        return M
    
    @staticmethod
    def AssembleForceVector( target_fun, basis, domain ):
        basis_dim = len( basis )
        F = numpy.zeros( basis_dim )
        for i in range( 0, basis_dim ):
            F[i] = scipy.integrate.quad( lambda x: basis[i](x) * target_fun(x), domain[0], domain[1], limit=1000, epsrel=1e-12 )[0]
        return F

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

## L2 Error
def ComputeL2Error( fun1, fun2, domain ):
    variate = list( fun1.atoms( sympy.Symbol ) )[0]
    integrand = sympy.lambdify( variate, ( fun1 - fun2 )**2.0 )
    integral = scipy.integrate.quad( integrand, domain[0], domain[1] )
    l2_error = numpy.sqrt( integral )
    return l2_error

def ComputeL2ErrorFast( fun1, fun2, domain ):
    integrand = lambda x: ( fun1(x) - fun2(x) )**2.0
    integral = scipy.integrate.quad( integrand, domain[0], domain[1], limit=1000, epsrel=1e-12 )
    l2_error = numpy.sqrt( integral[0] )
    return l2_error

## UTILITIES
def BasisToLatexString( basis_array, basis_id_string ):
    latex_string = basis_id_string + "=" + "\\begin{bmatrix}"
    for basis in basis_array:
        latex_string += f"{sympy.latex( sympy.simplify( basis.as_expr() ) )} \\\\"
    latex_string += "\\end{bmatrix}"
    return latex_string
