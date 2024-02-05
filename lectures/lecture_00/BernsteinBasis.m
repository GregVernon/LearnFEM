function bFun = BernsteinBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [0 1] );
    bFun = sym( zeros( degree + 1, 1 ) );
for a=0:degree
    bFun(a+1) = nchoosek( degree, a ) * ( variate ^ a ) * ( ( 1 - variate ) ^ ( degree - a ) );
end
bFun = simplify( bFun );
end