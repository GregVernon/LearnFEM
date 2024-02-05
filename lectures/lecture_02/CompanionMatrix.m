function C = CompanionMatrix( polyfun )
c = fliplr( coeffs( polyfun, 'all' ) );
C = sym( zeros( length(c) - 1 ) );
C(2:end, 1:end-1) = eye( length(c) - 2 );
C(:, end) = - ( c( 1:end-1 ) ./ c(end) );
end