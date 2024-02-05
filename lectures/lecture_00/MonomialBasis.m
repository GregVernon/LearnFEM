function bFun = MonomialBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [0, 1] );
bFun = sym( zeros( degree + 1, 1 ) );
for ii = 1 : degree + 1
    bFun(ii) = variate ^ (ii-1);
end
bFun = simplify( bFun );
end