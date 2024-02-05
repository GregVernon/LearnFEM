function basisFun = PolynomialBasisFunction( basisName, degree, variate, domain )
switch basisName
    case "Bernstein"
        basisFun = BernsteinBasis( degree, variate, domain );
    case "Chebyshev"
        basisFun = ChebyshevBasis( degree, variate, domain );
    case "Lagrange"
        basisFun = LagrangeBasis( degree, variate, domain );
    case "Lagrange-Chebyshev"
        basisFun = LagrangeChebyshevBasis( degree, variate, domain );
    case "Lagrange-Legendre"
        basisFun = LagrangeLegendreBasis( degree, variate, domain );
    case "Legendre"
        basisFun = LegendreBasis( degree, variate, domain );
    case "Monomial"
        basisFun = MonomialBasis( degree, variate, domain );
end
end