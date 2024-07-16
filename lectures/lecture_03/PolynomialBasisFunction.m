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
    case "Lagrange-NonUniform"
        basisFun = LagrangeNonUniformBasis( degree, variate, domain );
    case "Legendre"
        basisFun = LegendreBasis( degree, variate, domain );
    case "Monomial"
        basisFun = MonomialBasis( degree, variate, domain );
end
end

%% BASIS FUNCTIONS

function bFun = BernsteinBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [0 1] );
bFun = sym( zeros( degree + 1, 1 ) );
for a=0:degree
    bFun(a+1) = nchoosek( degree, a ) * ( variate ^ a ) * ( ( 1 - variate ) ^ ( degree - a ) );
end
end

function bFun = ChebyshevBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [-1, 1] );
bFun = sym( zeros( degree + 1, 1 ) );
for ii = 0 : degree
    if ii == 0
        bFun(ii+1) = variate ^ 0;
    elseif ii == 1
        bFun(ii+1) = variate ^ 1;
    else
        bFun(ii+1) = ( 2 * variate * bFun(ii) ) - bFun(ii-1);
    end
end
end

function bFun = MonomialBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [0, 1] );
bFun = sym( zeros( degree + 1, 1 ) );
for ii = 1 : degree + 1
    bFun(ii) = variate ^ (ii-1);
end
end

function bFun = LagrangeBasis( degree, variate, domain )
node = linspace( domain(1), domain(2), degree + 1 );
bFun = sym( zeros( degree+1, 1 ) );
for ii=1:degree+1  % ii is the current nodal basis function we're building
    bFun(ii) = variate ^ 0;
    for jj = 1 : degree + 1 % jj is evaluating the product series for the current node
        if ii ~= jj
            bFun(ii) = bFun(ii) * ( ( variate - node(jj) ) / ( node(ii) - node(jj) ) );
        end
    end
end
end

function bFun = LagrangeChebyshevBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [-1, 1] );
T = ChebyshevBasis( degree + 1, variate, [-1, 1] );
node = real( vpa( PolynomialRoots( T(end) ) ) );
bFun = sym( zeros( degree + 1, 1 ) );
for ii = 1 : degree + 1  % ii is the current nodal basis function we're building
    bFun(ii) = variate ^ 0;
    for jj = 1 : degree + 1 % jj is evaluating the product series for the current node
        if ii ~= jj
            bFun(ii) = bFun(ii) * ( (variate - node(jj) ) / ( node(ii) - node(jj) ) );
        end
    end
end
end

function bFun = LagrangeLegendreBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [-1, 1] );
L = LegendreBasis( degree + 1, variate, [-1, 1] );
node = real( vpa( PolynomialRoots( L(end) ) ) );
bFun = sym( zeros(degree+1,1) );
for ii = 1 : degree + 1  % ii is the current nodal basis function we're building
    bFun(ii) = variate ^ 0;
    for jj = 1 : degree + 1 % jj is evaluating the product series for the current node
        if ii ~= jj
            bFun(ii) = bFun(ii) * ( (variate - node(jj) ) / ( node(ii) - node(jj) ) );
        end
    end
end
end

function bFun = LegendreBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [-1, 1] );
bFun = sym( zeros( degree + 1, 1 ) );
for ii = 0 : degree
    if ii == 0
        bFun(ii+1) = variate ^ 0 ;
    elseif ii == 1
        bFun(ii+1) = variate ^ 1;
    else
        n = ii - 1;
        term1 = ( 2 * n + 1 ) * variate * bFun(n+1);
        term2 = n * bFun(n);
        bFun(ii+1) = simplify ( ( term1 - term2 ) / ( n + 1 ) );
    end
end
end