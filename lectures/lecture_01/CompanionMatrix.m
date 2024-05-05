function C = CompanionMatrix( polyfun )
c = fliplr( coeffs( polyfun, 'all' ) );
if length( c ) > 2
    C = sym( zeros( length(c) - 1 ), "r" );
    C(2:end, 1:end-1) = eye( length(c) - 2 );
    C(:, end) = - ( c( 1:end-1 ) ./ c(end) );
elseif length( c ) == 2
    C =  c(1) / (-c(2)) ;
elseif isscalar( c )
    C = [];
end