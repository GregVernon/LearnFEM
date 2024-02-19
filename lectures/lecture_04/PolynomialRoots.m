function [poly_roots, C] = PolynomialRoots( polyfun )
c = fliplr( coeffs( polyfun, 'all' ) );
C = sym( zeros( length(c) - 1 ) );
C(2:end, 1:end-1) = eye( length(c) - 2 );
C(:, end) = - ( c( 1:end-1 ) ./ c(end) );
poly_roots = eig( C );
[~, sidx] = sort( double ( real( vpa( poly_roots ) ) ) );
poly_roots = poly_roots( sidx );
end