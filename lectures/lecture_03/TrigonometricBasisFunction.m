function basisFun = TrigonometricBasisFunction( basisName, degree, variate, domain )
switch basisName
    case "Trigonometric"
        basisFun = TrigonometricBasis( degree, variate, domain );
    case "Exponential"
end
end

function bFun = TrigonometricBasis( degree, variate, domain )
variate = ChangeOfVariable( variate, domain, [0 pi] );
bFun = sym( zeros( 2 * degree + 1 , 1 ) );
bFun(1) = variate^0;
n = 1;
for a=1:degree
    n = n + 1;
    bFun(n) = sin( a * variate );
    n = n + 1;
    bFun(n) = cos( a * variate );
end
end