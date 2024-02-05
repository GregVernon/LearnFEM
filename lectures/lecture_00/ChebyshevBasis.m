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
bFun = simplify( bFun );
end