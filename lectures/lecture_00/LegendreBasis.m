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
bFun = simplify( bFun );
end