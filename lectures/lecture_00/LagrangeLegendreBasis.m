function bFun = LagrangeLegendreBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [-1, 1] );
L = LegendreBasis( degree + 1, variate, [-1, 1] );
node = simplify( real( PolynomialRoots( L(end) ) ) );
node = ChangeOfVariable( node, [min( node ), max( node )], [-1, 1] );
bFun = sym( zeros(degree+1,1) );
for ii = 1 : degree + 1  % ii is the current nodal basis function we're building
    bFun(ii) = variate ^ 0;
    for jj = 1 : degree + 1 % jj is evaluating the product series for the current node
        if ii ~= jj
            bFun(ii) = bFun(ii) * ( (variate - node(jj) ) / ( node(ii) - node(jj) ) );
        end
    end
end
bFun = simplify( expand( bFun ) );
end