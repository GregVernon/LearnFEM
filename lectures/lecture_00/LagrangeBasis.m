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
bFun = simplify( bFun );
end