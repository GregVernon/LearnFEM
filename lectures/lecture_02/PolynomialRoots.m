function poly_roots = PolynomialRoots( polyfun )
C = CompanionMatrix( polyfun );
poly_roots = eig( C );
[~, sidx] = sort( double ( real( vpa( poly_roots ) ) ) );
poly_roots = poly_roots( sidx );
end